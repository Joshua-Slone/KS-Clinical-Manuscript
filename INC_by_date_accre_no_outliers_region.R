inc_fxn <- function(dat,dat1,dat2,dat3,n){
  
  betas_gr_ea = vcovs_gr_ea = betas_ehr_ea = vcovs_ehr_ea = betas_gr_la = vcovs_gr_la = betas_ehr_la = vcovs_ehr_la = vector("list", n)
    
  for(i in 1:n){
    
    d6=dat
    d6_val=dat1
    d_1=dat2
    d1=dat3
    
    #d1 = d1 |> filter(!(patient %in% c("74605","pe.3397")))
    #d6 = d6 |> filter(!(patient %in% c("74605","pe.3397")))
    #d6_val = d6_val |> filter(!(patient %in% c("74605","pe.3397")))
    #d6_val = d6_val |> filter(!(patient %in% c("pe.3397")))
    d6_val = d6_val |> filter(!(patient %in% c("74605")))
    
    #d_1 = d_1 |> filter(!(patient %in% c("74605","pe.3397")))
    
    b1 = d6[d6$row==1,]
    ### imputation for phase 1 cd4 values based on error prone variables
    mod5 = lm(cd4_rt ~  male +program+ ns(age,df=4) + year +art +  ks*ns(time,df=2), data= b1) 
    ## to easily get design matrix for the full cohort
    d6_1 = b1 |> mutate(cd4_rt=1)
    v3=model.matrix(terms(mod5),d6_1)
    beta2 = mod5$coefficients
    vcov2 = vcov(mod5)
    rmvs2=rmvnorm(n=1,beta2,vcov2)
    rmvs2 = as.vector(rmvs2)
    cd4_rt_new = as.vector((v3 %*% rmvs2) + sample(resid(mod5),1))
    
    b1$cd4_rt_init_imp=if_else(is.na(b1$cd4_rt),cd4_rt_new,b1$cd4_rt)
    b2 = d6[d6$row!=1,] |> mutate(cd4_rt_init_imp=cd4_rt)
    d6=rbind(b1,b2) |> arrange(patient,row) |> mutate(cd4_rt_init_imp=if_else(cd4_rt_init_imp<=0,1,cd4_rt_init_imp)) |> fill(cd4_rt_init_imp,age)
    
    
    ### imputation for phase 2 cd4 values
    
    b1_incommon = b1 |> filter(patient %in% d6_val$patient) 
    
    b1_val = d6_val |> filter(row==1) |> left_join(b1_incommon[,c("patient","cd4_rt_init_imp")]) |> mutate(ch=abs(cd4_rt_init_imp-cd4_rt_final)) |> arrange(ch)
    
    mod6 = lm(cd4_rt_final ~  male_final +program+ ns(age_final,df=4) + ns(cd4_rt_init_imp,df=4) + year +art_final + ks_final*ns(time,df=2), data= b1_val)
    ## to easily get design matrix for phase 2
    d6_1 = b1_val |> mutate(cd4_rt_final=1)
    v4=model.matrix(terms(mod6),d6_1)
    beta2 = mod6$coefficients
    vcov2 = vcov(mod6)
    rmvs2=rmvnorm(n=1,beta2,vcov2)
    rmvs2 = as.vector(rmvs2)
    cd4_rt_final_new = as.vector((v4 %*% rmvs2) + sample(resid(mod6),1))
    
    b1_val$cd4_rt_imp_final=if_else(is.na(b1_val$cd4_rt_final),cd4_rt_final_new,b1_val$cd4_rt_final)
    b2_val = d6_val[d6_val$row!=1,] |> mutate(cd4_rt_imp_final=cd4_rt_final)
    b1_val = b1_val |> mutate(ch=NULL, cd4_rt_init_imp=NULL)
    d6_val=rbind(b1_val,b2_val) |> arrange(patient,row) |> mutate(cd4_rt_imp_final=if_else(cd4_rt_imp_final<=0,1,cd4_rt_imp_final)) |> fill(cd4_rt_imp_final)
    
    
    #rm(mod5, mod6)
    
    d61 = d6 |> mutate(age=round(as.numeric(age/365.25),3),fu=fu/365.25,cd4_rt_init_imp=round(cd4_rt_init_imp,3))
    
    #mod_ehr_inc = glm(event ~ male+program+I(age/10)+I(year-2010)+art + cd4_rt_init_imp,offset = log(fu/1000),data= d61, family = poisson) 
    
    
    ## Get these datasets to match for overlap phase-1 and phase-2, using old code or new code 
    d_c1=d6 |> filter(patient %in% unique(d6_val$patient)) |> mutate(min_fu=date,max_fu=as.Date(date+fu-1))
    d6_val = d6_val |> mutate(min_fu=date,max_fu=as.Date(date+fu-1))
    d6 = d6 |> mutate(min_fu=date,max_fu=as.Date(date+fu-1))
    
    
    ### Get overlapping timeframes
    ## get date ranges for each patient who was validated based on their phase 1 and phase 2 data
    ### phase 1
    cc1= d_c1 |> lazy_dt() |> group_by(patient) |> dplyr::summarize(n1=max(row_number()),d_min1=min(min_fu),d_max1=max(max_fu)) |> ungroup() |> as.data.frame()
    ### phase 2
    cc2= d6_val |> lazy_dt() |> group_by(patient) |> dplyr::summarize(n2=max(row_number()),d_min2=min(min_fu),d_max2=max(max_fu)) |> ungroup() |> as.data.frame()
    
    ### phase-1 and phase-2 date ranges and rows for V==1
    cc3=cc1 |> left_join(cc2) |> mutate(min_date1= pmax(d_min1,d_min2),max_date1= pmin(d_max1,d_max2))
    
    
    #################################################################################################################################
    #################################################################################################################################
    #################################################################################################################################
    ### Create overalapping phase-2 dataset for incidence analyses
    
    
    
    event_match <- cc3 |> select(patient,min_date1,max_date1) |>  pivot_longer(
      cols = c(min_date1,max_date1), 
      names_to = "date_type",
      values_to = "date") |> 
      left_join(cc3[,c('patient','min_date1','max_date1')])  
    
    merged_match <- event_match |>
      full_join(d6_val[,c("patient", "date","art_final", "cd4_rt_imp_final", "event_final")], by = c("patient", "date")) |>
      full_join(d_c1[,c("patient", "date","art", "cd4_rt_init_imp", "event")], by = c("patient", "date")) |>
      arrange(patient, date) |>
      group_by(patient) |> fill(min_date1,max_date1,.direction="updown") |> fill(art, art_final, cd4_rt_init_imp, cd4_rt_imp_final,event, event_final) |> filter(date >=min_date1 & date <= max_date1) |> dplyr::mutate(event_final=if_else(row_number() != n()-1,0,event_final),
                                                                                                                                                                                                                        start_date = min(date),start = as.numeric(date - start_date),
                                                                                                                                                                                                                        stop = lead(start)) |> filter(row_number() != n()) |>  mutate(year=year(date)) |> group_by(patient, year, art, art_final, cd4_rt_init_imp,cd4_rt_imp_final,event, event_final) |> 
      dplyr::summarize(min=min(start),max=max(stop),date=min(date),ct=n()) |>
      mutate(fu=max-min) |> 
      arrange(patient,date) |> ungroup() 
    
    
    d6_2_overlap = merged_match |> left_join(d1[,c('patient','male','male_final','program','birth_d','birth_d_final')]) |> group_by(patient) |> dplyr::mutate(age=if_else(row_number()==1,date-birth_d,NA),age_final=if_else(row_number()==1,date-birth_d_final,NA),fu=if_else(fu==0,1,fu)) |> fill(age,age_final) |> ungroup()
    

    d6_all= bind_rows(d6_2_overlap,d6[!(d6$patient %in% d_c1$patient),]) |> lazy_dt() |> group_by(patient) |> fill(age,age_final) |> mutate(V=if_else(patient %in% unique(d6_2_overlap$patient),1,0)) |> ungroup() |> as.data.frame()
    
    #################################################################################################################################
    #################################################################################################################################
    #################################################################################################################################
    
    ### EA
    
    d_2 = d_1 |> filter(patient %in% unique(d6_all$patient) & program==1)
    mod_w = glm(R ~ stop_d + program + death_y, family=binomial,data = d_2)
    mod_p<- mod_w$fitted.values
    d_p = data.frame(d_2,p_est=mod_p)
    
    #load(file="~/Desktop/PhD/Research/KS_Data/Code/wgts_by_hand.Rda")
    load("strata.Rda")
    
    #load(file="~/Desktop/PhD/Research/KS_Data/Code/strata.Rda")
    
    #################################################################################################################
    #### estimated weights including nonresponse probabilities from the incidence analysis cohort
    #load("/Users/joshslone/Library/CloudStorage/OneDrive-Vanderbilt/KS_MI_datasets/Base_data/d_p.Rda") 
    pb1=d6_all |> filter(program==1) |> left_join(d_c2,by="patient")|> dplyr::count(patient,strata,V) |> left_join(d_p[,c('patient','p_est')],by="patient")
    mod_w = glm(V ~ strata, family=binomial,data = pb1)
    mod_p <- mod_w$fitted.values
    pb2 = data.frame(pb1,p_est_overall=mod_p)
    pb2 = pb2 |> mutate(w_est = 1/(p_est_overall*p_est))
    #d_inc_cw = p_all |> left_join(pb2[,c('patient','w_est')],by="patient") |> mutate(fu=as.numeric(as.Date(max_fu)-as.Date(min_fu)+1))
    
    
    #d_inc_cw = d6_all |> left_join(pb2[,c('patient','w_est','strata')],by="patient") |> mutate(age=round(as.numeric(age/365.25),4),age_final=round(as.numeric(age_final/365.25),4),fu=round(fu/365.25,3),cd4_rt_init_imp=round(cd4_rt_init_imp,2),cd4_rt_imp_final=round(cd4_rt_imp_final,2))
    
    d_inc_cw = d6_all |> filter(program==1) |> left_join(pb2[,c('patient','w_est','strata')],by="patient") |> mutate(age=round(as.numeric(age/365.25),3),age_final=round(as.numeric(age_final/365.25),3),fu=fu/365.25,cd4_rt_init_imp=round(cd4_rt_init_imp,3),cd4_rt_imp_final=round(cd4_rt_imp_final,3))
    
    #################################################################################################################
    
    mod_ehr_inc_ea = glm(event ~ male+I(age/10)+I(year-2010)+art + cd4_rt_init_imp,offset = log(fu/1000),data= d_inc_cw, family = poisson)
    
    #### Raking
    IF_func <- function(fit){
      MM <- model.matrix(fit)
      Score.P   <- fit$y*MM - MM*exp(fit$linear.predictors)
      InfoMat.P <- t(MM*sqrt(exp(fit$linear.predictors))) %*% (MM*sqrt(exp(fit$linear.predictors)))/nrow(MM)
      inffunc <- (Score.P) %*% solve(InfoMat.P)
      inffunc
    }    
    naiveInfl3 =  IF_func(mod_ehr_inc_ea)[,-1]
    colnames(naiveInfl3)=paste("INF", 1:ncol(naiveInfl3), sep="")
    
    d14=cbind(d_inc_cw,naiveInfl3) 
    
    infs_eq1 <- as.formula(paste('~', paste0('INF', 1:ncol(naiveInfl3), collapse = ' + ')))
    designMF4 = twophase(id=list(~1,~1),weights=list(~1,~w_est),strata=list(~strata,NULL), data=d14,subset=~V==1,method="approx")
    Ninflcal3<-calibrate(designMF4,formula=infs_eq1,phase=2,calfun="raking") 
    fit_gr_inc_ea = svyglm(event_final ~ male_final + I(age_final/10) + I(year-2010) + art_final + cd4_rt_imp_final, offset = log(fu/1000), design=Ninflcal3, family = poisson)
    
    
    #################################################################################################################################
    #################################################################################################################################
    #################################################################################################################################
    
    ### LA
    
    d_2 = d_1 |> filter(patient %in% unique(d6_all$patient) & program==0)
    mod_w = glm(R ~ stop_d + program + death_y, family=binomial,data = d_2)
    mod_p<- mod_w$fitted.values
    d_p = data.frame(d_2,p_est=mod_p)
    
    #load(file="~/Desktop/PhD/Research/KS_Data/Code/wgts_by_hand.Rda")
    #load("strata.Rda")
    
    #################################################################################################################
    #### estimated weights including nonresponse probabilities from the incidence analysis cohort
    #load("/Users/joshslone/Library/CloudStorage/OneDrive-Vanderbilt/KS_MI_datasets/Base_data/d_p.Rda") 
    pb1=d6_all |> filter(program==0) |> left_join(d_c2,by="patient")|> dplyr::count(patient,strata,V) |> left_join(d_p[,c('patient','p_est')],by="patient")
    mod_w = glm(V ~ strata, family=binomial,data = pb1)
    mod_p <- mod_w$fitted.values
    pb2 = data.frame(pb1,p_est_overall=mod_p)
    pb2 = pb2 |> mutate(w_est = 1/(p_est_overall*p_est))
    #d_inc_cw = p_all |> left_join(pb2[,c('patient','w_est')],by="patient") |> mutate(fu=as.numeric(as.Date(max_fu)-as.Date(min_fu)+1))
    
    
    #d_inc_cw = d6_all |> left_join(pb2[,c('patient','w_est','strata')],by="patient") |> mutate(age=round(as.numeric(age/365.25),4),age_final=round(as.numeric(age_final/365.25),4),fu=round(fu/365.25,3),cd4_rt_init_imp=round(cd4_rt_init_imp,2),cd4_rt_imp_final=round(cd4_rt_imp_final,2))
    
    d_inc_cw = d6_all |> filter(program==0) |> left_join(pb2[,c('patient','w_est','strata')],by="patient") |> mutate(age=round(as.numeric(age/365.25),3),age_final=round(as.numeric(age_final/365.25),3),fu=fu/365.25,cd4_rt_init_imp=round(cd4_rt_init_imp,3),cd4_rt_imp_final=round(cd4_rt_imp_final,3))
    
    #################################################################################################################
    
    mod_ehr_inc_la = glm(event ~ male+I(age/10)+I(year-2010)+art + cd4_rt_init_imp,offset = log(fu/1000),data= d_inc_cw, family = poisson)
    
    naiveInfl3 =  IF_func(mod_ehr_inc_la)[,-1]
    colnames(naiveInfl3)=paste("INF", 1:ncol(naiveInfl3), sep="")
    
    d14=cbind(d_inc_cw,naiveInfl3) 
    
    #infs_eq1 <- as.formula(paste('~', paste0('INF', 1:ncol(naiveInfl3), collapse = ' + ')))
    designMF4 = twophase(id=list(~1,~1),weights=list(~1,~w_est),strata=list(~strata,NULL), data=d14,subset=~V==1,method="approx")
    Ninflcal3<-calibrate(designMF4,formula=~INF2 + INF3 + INF4 + INF5,phase=2,calfun="raking") 
    fit_gr_inc_la = svyglm(event_final ~ male_final + I(age_final/10) + I(year-2010) + art_final + cd4_rt_imp_final, offset = log(fu/1000), design=Ninflcal3, family = poisson)
    

    
    
    ###############################################################################################################################################
    ###############################################################################################################################################
    ### stored results for each loop run 
    
    betas_gr_ea[[i]]  = coef(fit_gr_inc_ea)
    vcovs_gr_ea[[i]]  = vcov(fit_gr_inc_ea)
    betas_ehr_ea[[i]] = coef(mod_ehr_inc_ea)
    vcovs_ehr_ea[[i]] = vcov(mod_ehr_inc_ea)
    
    betas_gr_la[[i]]  = coef(fit_gr_inc_la)
    vcovs_gr_la[[i]]  = vcov(fit_gr_inc_la)
    betas_ehr_la[[i]] = coef(mod_ehr_inc_la)
    vcovs_ehr_la[[i]] = vcov(mod_ehr_inc_la)
    
    
  }
  
  ddest <- list("coef_gr_ea" = betas_gr_ea,"vcovs_gr_ea" = vcovs_gr_ea,"coef_gr_la" = betas_gr_la,"vcovs_gr_la" = vcovs_gr_la,"coef_ehr_ea" = betas_ehr_ea,"vcovs_ehr_ea" = vcovs_ehr_ea,"coef_ehr_la" = betas_ehr_la,"vcovs_ehr_la" = vcovs_ehr_la)
  
  return(ddest) 
  
}






###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################





