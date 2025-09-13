

# Log Prior function ------------------------------------------------------


Logprior<-function(prior_dsl){
  
  function(theta){
    
    theta_exp<-exp(theta)
    
    monty::monty_model_density(prior_dsl, theta_exp)
  }
}



# Log Likelihood function -------------------------------------------------


my_Loglikelihood<-function(target_data,
                           p_list,
                           model,
                           update_p){

  
  function(theta){
    
    params<-exp(theta)
    
    
    # Set staring parameters
    pars<-update_p(params,p_list$pars_list)  
    
    # Determinsitic object
    sys  <- dust2::dust_system_create(model, pars, deterministic = TRUE, n_particles = 1)
    
    # Get state index
    index<- dust2::dust_unpack_index(sys)
    
    #Get default state initial
    state<-dust2::dust_system_state(sys)
    
    named_state<-dust2::dust_unpack_state(sys,state)
    
    # Set default initial state of the model 
    dust2::dust_system_set_state_initial(sys)
    
    # Set parameters
    dust2::dust_system_update_pars(sys, pars=pars)
    
    
    # run simulation
    endsim<-tail(target_data$time,1)
    tt<-seq(0,endsim)
    
    #tt<- seq(0,1000)
    
    
    
    y0<- dust2::dust_system_simulate(sys, tt)
    
    
    # unpack results
    state <- dust2::dust_unpack_state(sys, y0)

    
    
# 
#     matplot(t(state$S),type="l")
# 
# 
#     plot(colSums(state$pop_by4age),type="l")
# 
# 
#     plot((state$reported_wk_gi3),type="l")
# 
#     plot((state$reported_wk_gi),type="l")
# 
#     plot((state$reported_wk_gii4),type="l")
# 
#    plot((state$reported_wk_gii),type="l")
# 
#    # 
#    browser()
        
    exp_noise <- 1e6
    
    noise<-rexp(n = 1, rate = exp_noise)
    
    # IID2 data incidence -----------------------------------------------------
    
    model_t<-target_data$time[which(!is.na(target_data$cases_a1.1))]-1
    data_t <- which(!is.na(target_data$cases_a1.1))
    
    modelled_irate <-rbind(
      ( (  state$inc_year_gi3[1,model_t]+
             state$inc_year_gi[1,model_t]+
             state$inc_year_gii4[1,model_t]+
             state$inc_year_gii[1,model_t])/(state$pop_by4age[1,model_t])) + noise,
      
      ( (  state$inc_year_gi3[2,model_t]+
             state$inc_year_gi[2,model_t]+
             state$inc_year_gii4[2,model_t]+
             state$inc_year_gii[2,model_t])/(state$pop_by4age[2,model_t])) + noise,
      
      ( (  state$inc_year_gi3[3,model_t]+
             state$inc_year_gi[3,model_t]+
             state$inc_year_gii4[3,model_t]+
             state$inc_year_gii[3,model_t])/(state$pop_by4age[3,model_t])) + noise,
      
      ( (  state$inc_year_gi3[4,model_t]+
             state$inc_year_gi[4,model_t]+
             state$inc_year_gii4[4,model_t]+
             state$inc_year_gii[4,model_t])/(state$pop_by4age[4,model_t])) + noise,
      
      ( (  state$inc_year_gi3[5,model_t]+
             state$inc_year_gi[5,model_t]+
             state$inc_year_gii4[5,model_t]+
             state$inc_year_gii[5,model_t])/(state$pop_by4age[5,model_t])) + noise
    )*1000
    
    observed_size<-c(
      26.9,   # Person-years in 0 to 1
      190.8,  # Person-years in 1 to 4
      424.1,  # Person-years in 5 to 14
      2647.8, # Person-years in 15 to 65
      1369.1  # Person-years in 65plus
    )
    
    observations_irate <-rbind(
      target_data$cases_a1[data_t],
      target_data$cases_a1.1[data_t],
      target_data$cases_a1.2[data_t],
      target_data$cases_a1.3[data_t],
      target_data$cases_a1.4[data_t]
    )
    
    # llk_irate<-colSums(dbinom(x =  observations_irate,
    #                               size = round(target_data_size),
    #                               prob =  modelled_irate,
    #                               log = TRUE),na.rm=TRUE)
    
    
    llk_irate<-colSums(dpois(x =  round(observations_irate),
                             lambda = modelled_irate,
                             log = TRUE),na.rm=TRUE)
    
    
    
    # SGSS weekly reported series all -------------------------------------------------------
    
    model_t<-target_data$time[which(!is.na(target_data$reported))]-1
    data_t <- which(!is.na(target_data$reported))


    modelled_report <- colSums(state$reported_wk[,model_t]) +  noise

    observations_reported <- target_data$reported[data_t]

# 
#      negbin_dispersion<-params['reported_var']
#     # # #
#     llk_reported<-sum(dnbinom(x = observations_reported,
#                           size = modelled_report,prob= negbin_dispersion , log = TRUE))/100
#     
#     
    
    # =========Poisson
    # 
    llk_reported<-sum(dpois(x = observations_reported,
                        lambda= modelled_report, log = TRUE))/313


    # 
    # plot(observations_reported)
    # lines(modelled_report)

    # browser()
    #  modelled_report2 <- colSums(state$reported_wk[,])
    # 
    # plot(modelled_report2,xlim = c(1000,3000),ylim=c(0,1000),type="l")
    # 
    # 
    
    
    # SGSS weekly reported series average 5 seasons -------------------------------------------------------
    
    # model_t<-target_data$time[which(!is.na(target_data$reported_avg))]-1
    # data_t <- which(!is.na(target_data$reported_avg))
    # 
    # 
    # modelled_report <- colSums(state$reported_wk[,model_t]) +  noise
    # 
    # observations_reported <- target_data$reported_avg[data_t]
    # 
    # 
    # negbin_dispersion<-params['reported_var']
    # # # #
    # llk_reported<-sum(dnbinom(x = observations_reported,
    #                           size = modelled_report,prob= negbin_dispersion , log = TRUE))
    
    

    
    
    
    
    # SGSS weekly reported series by strain-------------------------------------------------------
    
    model_t<-target_data$time[which(!is.na(target_data$reported_gi3))]-1 # model strats from 0 and to match by week output
    data_t <- which(!is.na(target_data$reported_gi3))
    
    mod<-colSums(data.frame(
    modelled_report_gi3 = state$reported_wk_gi3[model_t]  + noise,
    modelled_report_gi  = state$reported_wk_gi[model_t]   + noise,
    modelled_report_gii4= state$reported_wk_gii4[model_t] + noise,
    modelled_report_gii = state$reported_wk_gii[model_t]  + noise))
    
    dat<-colSums(data.frame(
    observations_reported_gi3  = target_data$reported_gi3[data_t],
    observations_reported_gi   = target_data$reported_gi[data_t],
    observations_reported_gii4 = target_data$reported_gii4[data_t],
    observations_reported_gii  = target_data$reported_gii[data_t]))
    
    
    
    llk_reported_bystrain <- sum(dpois(x = dat,
                                  lambda= mod, log = TRUE))
    


    # browser()
    # totals<-colSums(mod)
    # 
    # props<-totals/sum(totals)
    # 
    # df<-data.frame(value=c(props),strain=names(props))
    # 
    # 
    # library(ggplot2)
    # ggplot(df,aes(x=strain,y=value, fill = strain))+
    #   geom_bar(stat = "identity")
    
    
    
    
    
    # negbin_dispersion<-pars$reported_var  
    # # # #  
    # llk_reported_gi3<-dnbinom(x = observations_reported_gi3,
    #                       size = modelled_report_gi3,prob= negbin_dispersion , log = TRUE)/52
    # 
    # llk_reported_gi<-dnbinom(x = observations_reported_gi,
    #                       size = modelled_report_gi,prob= negbin_dispersion , log = TRUE)/52
    # 
    # llk_reported_gii4<-dnbinom(x = observations_reported_gii4,
    #                       size = modelled_report_gii4,prob= negbin_dispersion , log = TRUE)/52
    # 
    # llk_reported_gii<-dnbinom(x = observations_reported_gii,
    #                       size = modelled_report_gii,prob= negbin_dispersion , log = TRUE)/52
    
    
    
    # llk_reported_gi3 <- sum(dpois(x = observations_reported_gi3,
    #                         lambda= modelled_report_gi3, log = TRUE))
    # 
    # llk_reported_gi  <- sum(dpois(x = observations_reported_gi,
    #                        lambda = modelled_report_gi, log = TRUE))
    # 
    # llk_reported_gii4<- sum(dpois(x = observations_reported_gii4,
    #                          lambda = modelled_report_gii4, log = TRUE))
    # 
    # llk_reported_gii <- sum(dpois(x = observations_reported_gii,
    #                         lambda = modelled_report_gii, log = TRUE))
    
    
    
    
    # Seroprev GII4 children 1 to 7 -------------------------------------------------------
    model_t<- target_data$time[which(!is.na(target_data$sero1))]-1 
    data_t <- which(!is.na(target_data$sero1))
    
    
    modelled_sero<-state$seroprev_gii4[2:7,model_t] + noise
    
    observed_event<-rbind(
      target_data$sero1[data_t],
      target_data$sero2[data_t],
      target_data$sero3[data_t],
      target_data$sero4[data_t],
      target_data$sero5[data_t],
      target_data$sero6[data_t]
    )
    
    observed_size<-rbind(
      103,
      107,
      121,
      124,
      122,
      109
    )
    
    llk_sero<-colSums(dbinom(x=round(observed_event),
                             size = observed_size,
                             prob = modelled_sero,
                             log = TRUE),na.rm = TRUE)
    
    
    
    # SGSS age reporting fractions -------------------------------------------------------
    model_t<- target_data$time[which(!is.na(target_data$a0_event))]-1 
    data_t <- which(!is.na(target_data$a0_event))
    
    denom <- sum(
      state$infections_day_gi3[,model_t]+
      state$infections_day_gi[,model_t]+
      state$infections_day_gii4[,model_t]+
      state$infections_day_gii[,model_t])
      
    
    
    
    # # Weekly cases reported 
    model_by_age<-(c(
      (state$infections_day_gi3[1,model_t]+
          state$infections_day_gi[1,model_t]+
          state$infections_day_gii4[1,model_t]+
          state$infections_day_gii[1,model_t]),
      
      (state$infections_day_gi3[2,model_t]+
         state$infections_day_gi[2,model_t]+
         state$infections_day_gii4[2,model_t]+
         state$infections_day_gii[2,model_t]),
      
      (state$infections_day_gi3[3,model_t]+
         state$infections_day_gi[3,model_t]+
         state$infections_day_gii4[3,model_t]+
         state$infections_day_gii[3,model_t]),
      
      (state$infections_day_gi3[4,model_t]+
         state$infections_day_gi[4,model_t]+
         state$infections_day_gii4[4,model_t]+
         state$infections_day_gii[4,model_t]))/ denom)+noise
      
 
   sample_size_fac<-12
 
   observed_size<-round(56308/sample_size_fac) 
    
   observed_event<-round(c(
     target_data$a0_event[data_t],
     target_data$a5_event[data_t],
     target_data$a15_event[data_t],
     target_data$a65_event[data_t]
   )) 
    
   llk_reported_by_age<- sum(dbinom(x =  observed_event,
                           size   =  observed_size,
                           prob   =  model_by_age,
                           log    = TRUE) )  
   
   
 
   
   # aggregated likelihood -------------------------------------------------------
   llk<-colSums(
     rbind(                  llk_irate,
                             llk_reported,
                             llk_reported_bystrain, 
                             # llk_reported_gi3,
                             # llk_reported_gi,
                             # llk_reported_gii4,
                             # llk_reported_gii,
                             llk_reported_by_age,
                             llk_sero
     ) , na.rm=T)


   return(llk) 
   
   
  }  
}