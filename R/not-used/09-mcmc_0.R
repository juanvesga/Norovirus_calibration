# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_dataS   <- file.path(root,"output", "data_short.qs2")
infile_dataL   <- file.path(root,"output", "data_long.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")
infile_prior   <- file.path(root,"output", "priors.qs2")
infile_limits  <- file.path(root, "output", "limits.qs2")
infile_prior2  <- file.path(root, "output", "priors2.qs2")
infile_limits2 <- file.path(root, "output", "limits2.qs2")
infile_model0  <- file.path(root,"models", "model_simple_fits.R")
infile_model1  <- file.path(root,"models", "model_simple_no_reinf.R")
infile_model2  <- file.path(root,"models", "model_2pars_fits.R")
infile_model3  <- file.path(root,"models", "model_drop_fits.R")
infile_index   <- file.path(root, "output", "index.qs")
infile_odin1chain   <- file.path(root, "output", "chains_imm2par.RData")


#outfile <- file.path(root, "output", "mcmc_chains.qs2")



#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))
source(file.path(root, "R", "likelihood_functions.R"))
# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(odin2)
library(dust2)
library(ggplot2)
library(dplyr)
library(fitR)
library(foreach)
library(doSNOW)


# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------

#model_takes= 0.0130696 #0minutes per run 
#model_takes<- 0.01414531 #1minutes per run 
# model_takes xxx #2minutes per run 
model_takes= 0.04073333 #3minutes per run 

## Cross protection 5%

model       <- "_0" # Simple SEIAR
#model       <- "_1" # Simple SEIAR no reinf no cross-prot
#model       <- "_2" # Full version with reinf and cross-protection
#model      <- "_3" # Full version with reinf and cross-protection and drop immunity


## Cross protection 25%

#model      <- "_4" # Simple SEIAR
#model      <- "_5" # Simple SEIAR no reinf no cross-prot
#model      <- "_6" # Full version with reinf and cross-protection
#model      <- "_7" # Full version with reinf and cross-protection and drop immunity


## Cross protection 50%

#model       <- "_8" # Simple SEIAR
#model      <- "_9" # Simple SEIAR no reinf no cross-prot
#model      <- "_10" # Full version with reinf and cross-protection
#model      <- "_11" # Full version with reinf and cross-protection and drop immunity

outfile <- file.path(root, "output", paste0("mcmc_chains",model,".qs2"))
outfile2 <- file.path(root, "output", paste0("cov",model,".qs2"))



# MCMC settings -----------------------------------------------------------
nchains    <- 4

niterations<- 40000  
use_last   <- TRUE # TRUE starts from previous run ans uses previous covmat 
size_adapt <- NULL # when to start proposal size adaptation (NULL is no adapt)
shape_start<- NULL # when to start proposal shape adaptation (NULL is no adapt)
paralell   <- 1    # 1 = run 4 chains in parallel; 0 = run 4 chains sequentially 

exp_t<-model_takes*niterations

print(c(exp_t," mins and",exp_t/60, "hours" ))

# Source Odin model
if (model=="_0"||model=="_4"||model=="_8"){
  
  source(infile_model0)
  
  prior      <- qs_read(infile_prior2)
  limits     <- qs_read(infile_limits2)
  
  if (model=="_0"){
    crossp <-0.05
  }else if (model=="_4"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
}else if((model=="_1"||model=="_5"||model=="_9")){
  
  source(infile_model1)
  
  prior      <- qs_read(infile_prior2)
  limits     <- qs_read(infile_limits2)
  
  if (model=="_1"){
    crossp <-0.05
  }else if (model=="_5"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
}else if((model=="_2"||model=="_6"||model=="_10")){
  source(infile_model2)
  
  prior      <- qs_read(infile_prior)
  limits     <- qs_read(infile_limits)
  
  if (model=="_2"){
    crossp <-0.05
  }else if (model=="_6"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
}else if((model=="_3"||model=="_7"||model=="_11")){
  source(infile_model3)
  
  prior      <- qs_read(infile_prior2)
  limits     <- qs_read(infile_limits2)
  
  if (model=="_3"){
    crossp <-0.05
  }else if (model=="_7"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
  
}


#parameters <- qs_read(infile0)



data       <- qs_read(infile_dataS)
pars_list  <- qs_read(infile_input)
pars_list$pars_list$crossp_GI=crossp
pars_list$pars_list$crossp_GII=crossp

#Create LogPrior handle
my_prior<-Logprior(prior)



#Create logllk handle
compare_model<-my_Loglikelihood(data,
                                pars_list,
                                noro_model(),
                                update_parameters)




posterior<-function(fn_llk,
                    fn_prior){
  function(theta){
    fn_prior(theta)+fn_llk(theta)
  }
}

my_posterior<-posterior(compare_model,my_prior)




if (use_last){
  #combine chains
  chains     <- qs_read(outfile) 
  covs       <- qs_read(outfile2)
  longchain<-do.call(rbind, chains)
  
  
  #get last sample
  params<- rbind(
    unlist(exp(chains[[1]][nrow(chains[[1]]),1:(ncol(chains[[1]])-1)])),
    unlist(exp(chains[[2]][nrow(chains[[2]]),1:(ncol(chains[[2]])-1)])),
    unlist(exp(chains[[3]][nrow(chains[[3]]),1:(ncol(chains[[3]])-1)])),
    unlist(exp(chains[[4]][nrow(chains[[4]]),1:(ncol(chains[[4]])-1)])))
  
  
  # Build a dense maytyrix using empirical covmat
  
  
  epsilon = 1e-6
  
  d = nrow(covs[[1]])
  
  
  # 3. Create the stability matrix: epsilon * I
  stability_matrix = epsilon * diag(d)
  
  covmatrix<- list(
    (2.38^2 / d) * covs[[1]] +   stability_matrix,
    (2.38^2 / d) * covs[[2]] +   stability_matrix,
    (2.38^2 / d) * covs[[3]] +   stability_matrix,
    (2.38^2 / d) * covs[[4]] +   stability_matrix)
  
  
  # covmatrix<- list(
  #    covs[[1]],
  #    covs[[2]],
  #    covs[[3]],
  #    covs[[4]])
  
  
  
  
  
} else {
  
  if((model=="_0"||model=="_1"||model=="_4"||model=="_5"||model=="_8"||model=="_9")){
    
    # start of model immunity drop                                
    tmp<-c(
      beta_1      =    9.441207e-02,   
      beta_2      =    9.425697e-02 , 
      beta_3      =    1.899702e-01,  
      beta_4      =    9.409377e-02,  
      aduRR       =    1.051449e-01,  
      maternalAB  =    9.895775e+01,  
      imm_yr      =    6.420126e+01, 
      repfac_0    =     2.379536e+02, 
      repfac_5    =    1.412386e+03,  
      repfac_15   =    6.140188e+02,  
      repfac_65p  =    4.101715e+01,   
      #  reported_var=    4.863380e-01,  
      season_lag  =    1.907634e+01,  
      season_amp  =    1.052722e-01 )
    
  }else if((model=="_2"||model=="_6"||model=="_10")){
    tmp<-c(
      beta_1      =    0.09933151, 
      beta_2      =    0.09937132, 
      beta_3      =    0.2537905, 
      beta_4      =    0.1117861, 
      aduRR       =    0.1137812,    
      maternalAB  =    94.4391, 
      imm_yr      =    58.77856, 
      imm_fac     =    2.025501, 
      repfac_0    =    233.0642, 
      repfac_5    =    1319.582,  
      repfac_15   =    561.4284, 
      repfac_65p  =    32.78976,   
      #reported_var=    4.889772e-01, 
      season_lag  =    11.12073,  
      season_amp  =    0.1140431 )
    
  } else if((model=="_3"||model=="_7"||model=="_11")){
    
    # start of model immunity drop                                
    tmp<-c(
      beta_1      =    8.417229e-02, 
      beta_2      =    8.370465e-02, 
      beta_3      =    2.216278e-01, 
      beta_4      =    9.375149e-02, 
      aduRR       =    2.480987e-01, 
      maternalAB  =    1.553466e+02, 
      imm_yr      =    3.212510e+01,
      repfac_0    =    2.303229e+02,
      repfac_5    =    1.454593e+03, 
      repfac_15   =    6.573526e+02, 
      repfac_65p  =    3.952568e+01, 
      # reported_var=    4.889772e-01, 
      season_lag  =    1.110333e+01, 
      season_amp  =    1.052237e-01)
  }
  
  
  
  
  params<-rbind(tmp,tmp,tmp,tmp)
  rownames(params)<-NULL
  
  # get cov matrix from previous run
  covs<- abs(diag(log(tmp)))*0.000001#cov(chains_odin)
  colnames(covs)<-paste(names(tmp))
  rownames(covs)<-paste(names(tmp))
  
  covmatrix<-list(
    covs,
    covs,
    covs,
    covs)
  
  
}

# test it
my_prior(log(params[1,]))


compare_model(log(params[1,]))

my_posterior(log(params[1,]))




# MCMC foreach function ---------------------------------------------------



pb <- txtProgressBar(max = nchains, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foreach_fun<-function(nchains,
                      niterations,
                      params,
                      covmat,
                      limits,
                      posterior_fn,
                      size_adapt,
                      shape_start,
                      para,
                      .combine = 'c'){
  
  if (para==0){
    foreach(i = 1:nchains,.options.snow = opts) %do% {
      
      #size_adapt=NULL
      #shape_start=NULL
      
      res_mcmc <- fitR::mcmcMh(
        target=posterior_fn,
        initTheta=log(params[i,]),
        proposalSd = NULL,#diag(covmat),#log(params)/20,
        nIterations=niterations,
        covmat = covmat[[i]],
        limits = limits,
        adaptSizeStart = size_adapt,# niterations+400,
        adaptSizeCooling = 0.99,
        adaptShapeStart = shape_start,# niterations+500,
        adaptShapeStop = NULL,
        printInfoEvery = round(niterations*0.1),
        verbose = FALSE,
        maxScalingSd = 50
      ) 
      
      vc <- res_mcmc$trace
      trace <- apply(vc, 2, unlist)
      trace <- data.frame(trace)
      
      out<-list(
        trace=trace,
        acc_rate=res_mcmc$acceptanceRate,
        cov_mat =res_mcmc$covmatEmpirical)
    }
  }else{
    foreach(i = 1:nchains,.options.snow = opts) %dopar% {
      
      #size_adapt=NULL
      #shape_start=NULL
      
      res_mcmc <- fitR::mcmcMh(
        target=posterior_fn,
        initTheta=log(params[i,]),
        proposalSd = NULL,#diag(covmat),#log(params)/20,
        nIterations=niterations,
        covmat = covmat[[i]],
        limits = limits,
        adaptSizeStart = size_adapt,# niterations+400,
        adaptSizeCooling = 0.99,
        adaptShapeStart = shape_start,# niterations+500,
        adaptShapeStop = NULL,
        printInfoEvery = round(niterations*0.1),
        verbose = FALSE,
        maxScalingSd = 50
      ) 
      
      vc <- res_mcmc$trace
      trace <- apply(vc, 2, unlist)
      trace <- data.frame(trace)
      
      out<-list(
        trace=trace,
        acc_rate=res_mcmc$acceptanceRate,
        cov_mat =res_mcmc$covmatEmpirical)
    }
  }
  
  
  
  
}

# Call cores and register clusters

# ncores<-parallel::detectCores()
# 
# cl <- parallel::makeCluster(3)
# doParallel::registerDoParallel(cl)

cl <- makeCluster(nchains)
registerDoSNOW(cl)
start_time <- Sys.time()
results<-foreach_fun(nchains,
                     niterations,
                     params,
                     covmatrix,
                     limits,
                     my_posterior,
                     size_adapt,
                     shape_start,
                     paralell)

# Stop cluster
close(pb)
stopCluster(cl)


end_time <- Sys.time()
timelag<-end_time-start_time
print(timelag)
print(as.numeric(timelag)/niterations)


# post-processing ---------------------------------------------------------




traces<-rbind(results[[1]]$trace$logDensity,
              results[[2]]$trace$logDensity,
              results[[3]]$trace$logDensity,
              results[[4]]$trace$logDensity)

matplot(t(traces),type = "l")


# Coda
trace1 <- coda::mcmc(results[[1]]$trace)
trace2 <- coda::mcmc(results[[2]]$trace)
trace3 <- coda::mcmc(results[[3]]$trace)
trace4 <- coda::mcmc(results[[4]]$trace)

# combine traces as mcmc.list object
trace <- coda::mcmc.list(list(trace1, trace2,trace3,trace4))

library("lattice")  ## for the 'xyplot' command
xyplot(trace)



empiric_cov<-list(
  chain1=results[[1]]$cov_mat,
  chain2=results[[2]]$cov_mat,
  chain3=results[[3]]$cov_mat,
  chain4=results[[4]]$cov_mat
)


# Save results ------------------------------------------------------------



qs_save(trace,outfile)

qs_save(empiric_cov,outfile2)

