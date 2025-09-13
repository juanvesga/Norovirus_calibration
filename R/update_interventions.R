# Process model output function -------------------------------------------
update_parameters<-function(ins,params){

  
  if("imm_fac"%in%names(ins)){
  
  params$beta_1     = ins["beta_1"]
  params$beta_2     = ins["beta_2"]
  params$beta_3     = ins["beta_3"]
  params$beta_4     = ins["beta_4"]
  params$aduRR      = ins["aduRR"]
  params$maternalAB = ins["maternalAB"]
  params$imm_yr     = ins["imm_yr"]
  params$imm_fac    = ins["imm_fac"]
  params$repfac_0   = ins["repfac_0"]
  params$repfac_5   = ins["repfac_5"]
  params$repfac_15  = ins["repfac_15"]
  params$repfac_65p = ins["repfac_65p"]
  params$season_lag = ins["season_lag"]
  params$season_amp = ins["season_amp"]
  
  }else{
    
    params$beta_1     = ins["beta_1"]
    params$beta_2     = ins["beta_2"]
    params$beta_3     = ins["beta_3"]
    params$beta_4     = ins["beta_4"]
    params$aduRR      = ins["aduRR"]
    params$maternalAB = ins["maternalAB"]
    params$imm_yr     = ins["imm_yr"]
    params$repfac_0   = ins["repfac_0"]
    params$repfac_5   = ins["repfac_5"]
    params$repfac_15  = ins["repfac_15"]
    params$repfac_65p = ins["repfac_65p"]
    params$season_lag = ins["season_lag"]
    params$season_amp = ins["season_amp"]

    
  }
  
  
  
  return(params)
}


update_interventions<-function(mod,start_t,state0,run_t,param,label){
  
  # Set time   
  dust2::dust_system_set_time(mod, start_t)
  
  # Set model state  
  dust2::dust_system_set_state(mod,state0)
  
  # Set parameters
  dust2::dust_system_update_pars(mod, pars=param)
  
  # run simulation
  y0<- dust2::dust_system_simulate(mod, seq(start_t,start_t+run_t))
  
  # unpack results
  y0 <- dust2::dust_unpack_state(mod, y0)
  
  
  

  # Create output dataframe

  out<-data.frame(
    days                   = seq(0,run_t),  
    new_yearly_hosp        = y0$new_hosp_elder_yr+y0$new_hosp_adult_yr,
    new_weekly_cases       = y0$new_cases_week,
    new_weekly_deaths     =  y0$death_reported_wk,
    new_daily_cases        = y0$new_cases,
    cum_new_daily_cases    = cumsum(y0$new_cases),
    cum_new_weekly_cases   = cumsum(y0$new_cases_week),
    vacc_doses             = y0$n_vacc_doses,
    cum_vacc_doses         = cumsum(y0$n_vacc_doses),
    new_weekly_cases_gi3    = y0$new_cases_week_gi3,
    new_weekly_cases_gi     = y0$new_cases_week_gi,
    new_weekly_cases_gii4   = y0$new_cases_week_gii4,
    new_weekly_cases_gii    = y0$new_cases_week_gii,
    new_weekly_reported     = colSums(y0$reported_wk)
    
  )
  
  out$scenario<-label
  
  return(out)
  
}