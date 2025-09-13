# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_dataS   <- file.path(root,"output", "data_short.qs2")
infile_dataL   <- file.path(root,"output", "data_long.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")
infile_prior   <- file.path(root,"output", "priors.qs2")
infile_model2  <- file.path(root,"models", "model_2pars.R")
infile_chains  <- file.path(root, "output", "chains.qs")
infile_index   <- file.path(root, "output", "index.qs")
infile_odin1chain   <- file.path(root, "output", "chains_imm2par.RData")
outfile <- file.path(root, "output", "scenario_results.qs2")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(odin2)
library(dust2)
library(ggplot2)
library(dplyr)
library(foreach)

# Number of samples -------------------------------------------------------------

nsamps<- 500

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------


# Source Odin model
source(infile_model2)
parameters <- qs_read(infile0)
data       <- qs_read(infile_dataS)
pars_list  <- qs_read(infile_input)
prior      <- qs_read(infile_prior)
chains     <- qs_read(infile_chains) 
load(infile_odin1chain)

calibrated_pars<-exp(chains$pars)
#calibrated_pars<-exp(processed_chains$pars)
#init_states    <-processed_chains$state
#init_states    <-init_states[-c(120),]  ## Remove time vector from initial state

rm(chains)
rm(processed_chains)


# Sample from posterior
params<-calibrated_pars[sample(nrow(calibrated_pars), nsamps), ]#colMeans(calibrated_pars)
rm(calibrated_pars)


foreach_fun<-function(params,
                      pars_list,
                      noro_model,
                      update_parameters,
                      update_interventions,
                      .combine = 'c'){
  
  foreach(i = 1:nrow(params)) %dopar% {
    
    # Set staring parameters
    pars<-update_parameters(params[i,],pars_list$pars_list)  
    
    #dust object
    #sys  <- dust_system_create(noro_model(), pars, n_particles = 1)

    # Determinsitic object
    sys  <- dust2::dust_system_create(noro_model, pars, deterministic = TRUE, n_particles = 1)
    
    # Get state index
    index<- dust2::dust_unpack_index(sys)
    
    #Get default state initial
    state<-dust2::dust_system_state(sys)
    
    named_state<-dust2::dust_unpack_state(sys,state)
    
    # Set inital model state based on calibrated model
    
    # Set default initial state of the model 
    dust2::dust_system_set_state_initial(sys)
    
    # Set parameters
    dust2::dust_system_update_pars(sys, pars=pars)
    
    
    # run simulation
    endsim<-24109
    tt<-seq(0,endsim)
    y0<- dust2::dust_system_simulate(sys, tt)
    
    # unpack results
    y0 <- dust2::dust_unpack_state(sys, y0)
    
    
    # Get model state at the latest point for future simulations
    init_state0<-dust2::dust_system_state(sys)
    
    names(init_state0)<-paste0(names(y0))
    
    
    
    # Interventions -----------------------------------------------------------
    # index age
    #   1   0-1
    #   2   1-2
    #   3   2-3
    #   4   3-4
    #   5   5-6
    #   6   6-7
    #   7   7-15
    #   8   15-25
    #   9   25-35
    #   10  35-45
    #   11  45-55
    #   12  55-65
    #   13  65-75
    #   14  75-100
    
    # Scenario parameters
    vacc_eff         <- 0.75
    vacc_coverage    <- 0.9
    campaign_coverage<- 0.7
    t_zero           <- 1
    t_span           <-365*10
    t_campaign_stop  <-t_zero+180
    vacc_dis         <- 1
    vacc_trans       <- 0
    
    scen_tags<-c("Baseline",
                 "Routine under 1yr",
                 "Routine under 1 & 65yr",
                 "Routine under 1 + campaign under 5yr")
    
    
    # Baseline
    
    pars0<-list()
    
    runs0<-update_interventions(sys,t_zero,init_state0,t_span,pars0,"Baseline")
    runs0$cases_averted<-runs0$cum_new_daily_cases-runs0$cum_new_daily_cases
    runs0$dose_per_cases_averted<-runs0$cum_vacc_doses/(runs0$cum_new_daily_cases-runs0$cum_new_daily_cases)
    runs0$scenario<-scen_tags[1]
    
    # deaths<- as.data.frame(runs0$new_weekly_deaths) %>% 
    #   slice(seq(0, n(), by = 7))
    # 
    # 
    # admissions<- as.data.frame(runs0$new_yearly_hosp) %>% 
    #   slice(seq(0, n(), by = 365))
    # 
    # 
    # reported<-as.data.frame(runs0$new_weekly_reported) %>% 
    #   slice(seq(0, n(), by = 7))
    
    
    # Routine Vaccination for under 1
    
    pars1<-pars
    pars1$t_zero                 <- t_zero
    pars1$t_campaign_stop        <- t_campaign_stop 
    pars1$vaccination_coverage[1]<- vacc_eff*vacc_coverage
    pars1$vacc_switch_on[1]       <- 1
    pars1$campaign_switch[1]      <- 0
    pars1$vacc_dis                <- vacc_dis
    pars1$vacc_trans              <- vacc_trans
    
    
    runs1<-update_interventions(sys,t_zero,init_state0,t_span,pars1,"Routine under 1yr")
    runs1$cases_averted<-runs0$cum_new_daily_cases-runs1$cum_new_daily_cases
    runs1$dose_per_cases_averted<-runs1$cum_vacc_doses/(runs0$cum_new_daily_cases-runs1$cum_new_daily_cases)
    runs1$scenario<-scen_tags[2]
    
    # Routine Vaccination for under 1 and 65
    pars2<-pars
    pars2$t_zero                         <- t_zero
    pars2$t_campaign_stop                <- t_campaign_stop 
    pars2$vaccination_coverage[c(1,13)]  <- vacc_eff*vacc_coverage
    pars2$vacc_switch_on[c(1,13)]        <- 1
    pars2$campaign_switch[c(1,13)]       <- 0
    pars2$vacc_dis                       <- vacc_dis
    pars2$vacc_trans                     <- vacc_trans
    
    runs2<-update_interventions(sys,t_zero,init_state0,t_span,pars2,"Routine under 1 & 65yr")
    runs2$cases_averted<-runs0$cum_new_daily_cases-runs2$cum_new_daily_cases
    runs2$dose_per_cases_averted<-runs2$cum_vacc_doses/(runs0$cum_new_daily_cases-runs2$cum_new_daily_cases)
    runs2$scenario<-scen_tags[3]
    
    
    # Routine Vaccination for under 1 + 6mo campign under 5
    pars3<-pars
    pars3$t_zero                         <- t_zero
    pars3$t_campaign_stop                <- t_campaign_stop 
    pars3$vaccination_coverage[c(1:5)]   <- vacc_eff*vacc_coverage
    pars3$vacc_switch_on[c(1:5)]         <- 1
    pars3$campaign_switch[1]             <- 0
    pars3$campaign_switch[c(2:5)]        <- 1
    pars3$vacc_dis                       <- vacc_dis
    pars3$vacc_trans                     <- vacc_trans
    
    runs3<-update_interventions(sys,t_zero,init_state0,t_span,pars3,"Routine under 1 + campaign under 5yr")
    runs3$cases_averted<-runs0$cum_new_daily_cases-runs2$cum_new_daily_cases
    runs3$dose_per_cases_averted<-runs3$cum_vacc_doses/(runs0$cum_new_daily_cases-runs3$cum_new_daily_cases)
    runs3$scenario<-scen_tags[4]
    
    
    # Join together all results 
    ############# 

    
    
    df<-rbind(runs0,runs1,runs2,runs3)
    

    #df$scenario<- factor(df$scenario, levels=scen_tags)
    
  }
}

# Call cores and register clusters

ncores<-parallel::detectCores()

cl <- parallel::makeCluster(ncores-1)
doParallel::registerDoParallel(cl)

results<-foreach_fun(params,
                     pars_list,
                     noro_model(),
                     update_parameters,
                     update_interventions)

# Stop cluster
parallel::stopCluster(cl)



# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(results, outfile)


