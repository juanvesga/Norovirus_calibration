# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")
infile_model0  <- file.path(root,"models", "model_simple.R")
infile_model2  <- file.path(root,"models", "model_2pars.R")
infile_model3  <- file.path(root,"models", "model_drop.R")




#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))
source(file.path(root, "R", "collect_function.R"))
# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(abind, include.only = c("abind"))

library(odin2)
library(dust2)
library(ggplot2)
library(dplyr)
library(foreach)
library("lattice")  ## for the 'xyplot' command
library(fitR)


# Select model ------------------------------------------------------------


# Cross protection 5%

#model       <- "_0" # Simple SEIAR
#model       <- "_1" # Simple SEIAR no reinf no cross-prot
#model       <- "_2" # Full version with reinf and cross-protection
model      <- "_3" # Full version with reinf and cross-protection and drop immunity


# Cross protection 25%

#model      <- "_4" # Simple SEIAR
#model      <- "_5" # Simple SEIAR no reinf no cross-prot
#model      <- "_6" # Full version with reinf and cross-protection
#model      <- "_7" # Full version with reinf and cross-protection and drop immunity


# Cross protection 50%

#model       <- "_8" # Simple SEIAR
#model      <- "_9" # Simple SEIAR no reinf no cross-prot
#model      <- "_10" # Full version with reinf and cross-protection
#model      <- "_11" # Full version with reinf and cross-protection and drop immunity


# Number of samples -------------------------------------------------------------

nsamps<- 400
burn_in    <- 0
thin_in    <- 30



# Source Odin model
if (model=="_0"||model=="_4"||model=="_8"){
  
  source(infile_model0)

  
  if (model=="_0"){
    crossp <-0.05
  }else if (model=="_4"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
}else if((model=="_1"||model=="_5"||model=="_9")){
  
  source(infile_model1)
  

  
  if (model=="_1"){
    crossp <-0.05
  }else if (model=="_5"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
}else if((model=="_2"||model=="_6"||model=="_10")){
  source(infile_model2)

  
  if (model=="_2"){
    crossp <-0.05
  }else if (model=="_6"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
}else if((model=="_3"||model=="_7"||model=="_11")){
  source(infile_model3)
  

  
  if (model=="_3"){
    crossp <-0.05
  }else if (model=="_7"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
  
}

infile_chains <- file.path(root, "output", paste0("mcmc_chains",model,".qs2"))
outfile1      <- file.path(root, "output", paste0("mcmc_fits",model,".qs2"))
outfile2      <- file.path(root, "output", paste0("mcmc_inits",model,".qs2"))




# Source Odin model

pars_list  <- qs_read(infile_input)
pars_list$pars_list$crossp_GI=crossp
pars_list$pars_list$crossp_GII=crossp
trace     <- qs_read(infile_chains) 


# burn and thin
traceBurn <- burnAndThin(trace, burn = burn_in)
chains <- burnAndThin(traceBurn, thin = thin_in)


#combine chains

longchain<-do.call(rbind, chains)

id<-which(longchain[,"logDensity"]==max(longchain[,"logDensity"]))[1]

exp(longchain[id,])

#get last sample
params<- exp(longchain[sample(nrow(longchain),size=nsamps,replace=FALSE),])

# Functions ---------------------------------------------------------------

### Function to convert a list of dataframes to a 3D array

## All objects in the list will be dataframes with identical column headings.
list2ary = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}




foreach_fun<-function(params,
                      pars_list,
                      noro_model,
                      update_parameters,
                      collect_fn,
                      .combine='comb', 
                      .multicombine=TRUE){
  
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
    
    out<-collect_fn(y0)
    
    # # Get model state at the latest point for future simulations
    init_state0<-dust2::dust_system_state(sys)

    names(init_state0)<-paste0(names(y0))
    # 
    list(out,init_state0)
    
  }
}


# Execute runs ------------------------------------------------------------



# Call cores and register clusters

ncores<-parallel::detectCores()

cl <- parallel::makeCluster(ncores-1)
doParallel::registerDoParallel(cl)


results<-foreach_fun(params,
                     pars_list,
                     noro_model(),
                     update_parameters,
                     collect_function)

# Stop cluster
parallel::stopCluster(cl)



runs <- do.call(abind,c(lapply(results,function(x){x[[1]]}),along=3) )

runs <- aperm(runs, c(2,3,1))

initState <- as.data.frame(do.call(rbind,lapply(results,function(x){x[[2]]})))



qs_save(runs,outfile1)
qs_save(initState,outfile2)


