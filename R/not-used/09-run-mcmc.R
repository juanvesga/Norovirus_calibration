# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile0  <- file.path(root,"output", "parameters.qs2")
infile_dataS  <- file.path(root,"output", "data_short.qs2")
infile_dataL  <- file.path(root,"output", "data_long.qs2")
infile_input  <- file.path(root,"output", "params_list.qs2")
infile_prior  <- file.path(root,"output", "priors.qs2")
infile_model2  <- file.path(root,"models", "model_2pars.R")

outfile <- file.path(root, "output", "chains.qs")
outfile_index <- file.path(root, "output", "index.qs")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))


# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(odin2)
library(dust2)
library(monty)
library(posterior)
library(bayesplot)

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------
# MCMC conditions
nsamples=1000
nout    = 200 
nchains=3 
to_burn=nsamples*0.6
each_n = round((nsamples-to_burn)/nout) 
set.seed(42)

# Source Odin model
source(infile_model2)
parameters <- qs_read(infile0)
data       <- qs_read(infile_dataS)
pars_list  <- qs_read(infile_input)
prior      <- qs_read(infile_prior)
#Parameters object
pars<-pars_list$pars_list
packer<-pars_list$packer
start_list<-pars_list$start_pars



#dust object
sys  <- dust_system_create(noro_model(), pars, n_particles = 2)
index<- dust_unpack_index(sys)
state<-dust_system_state(sys)
named_state<-dust_unpack_state(sys,state)

# Generate Initial conditions vector
seed<-1# Number infected per strain
init<-generate_initial_state(index,state,parameters,seed,2)


# check initial conditions in current object
dust_system_state(sys)

# Pass initial conditions to model
dust_system_set_state_initial(sys)

# Check again initial conditions in object
dust_system_state(sys)

#Simulate system up to a point 
# t <- seq(0, 20000)
# y <- dust_system_simulate(sys, t)
# y <- dust_unpack_state(sys, y)
# 
# plot(t, y$new_cases_week, type = "l", xlab = "Time", ylab = "Infected population",
#      xlim = c(1e4,2e4),ylim=c(0,5e4))
# 
# plot(t[t %% 7 == 0], y$new_cases_week[t %% 7 == 0], type = "o", pch = 19,
#      ylab = "Infection incidence", xlab = "Time",xlim = c(1e4,10600),ylim=c(0,15e4))
# lines(t, y$new_cases_week, col = "red")



#dust_likelihood_run(filter_det, pars)

#stochastic
# filter     <- dust_filter_create(noro_model(), time_start = 0, data = data, n_particles=12)
# dust_likelihood_run(filter, pars)


# 
# likelihood <- dust_likelihood_monty(filter,
#                                     packer,
#                                     save_trajectories = TRUE)
# vcv <- diag(length(pars_list$start_pars)) * 0.00004
# sampler <- monty_sampler_random_walk(vcv)
# 
#  posterior <- likelihood + prior
# 
#  runner<-monty_runner_callr(4)
#  
#  samples <- monty_sample(posterior, 
#                          sampler, 
#                          50,
#                          initial = packer$pack(start_list), 
#                          n_chains = 1)
#  
 # 

# 
# # <- monty::monty_sampler_random_walk(vcv)
# 
# #vcv<-cov(samples_df[,1:14])
# sampler <- monty::monty_sampler_adaptive(vcv)
# 
# #runner <- monty::monty_runner_parallel(2)
# 
# # <- monty::monty_sample(posterior, sampler, 100,
#                       #   runner=runner,n_chains = 2)
# 
# samples <- monty::monty_sample(posterior, sampler, 200 ,n_chains = 1, initial = packer$pack(start_list))
# 
# 
# qs_save(samples, outfile)
# matplot(samples$density, type = "l", lty = 1,
#         xlab = "Sample", ylab = "Log posterior probability density")
# 
# 
# samples_df <- posterior::as_draws_df(samples)
# summary<-posterior::summarise_draws(samples_df)
# c(summary[,"median"])


# 
# bayesplot::mcmc_trace(samples_df)



### Improve covmatrix
# Create deterministic filter
 filter_det <- dust_unfilter_create(noro_model(), 
                                    time_start = 0, 
                                    data = data)
 
# 
#  likelihood <- dust_likelihood_monty(filter_det,
#                                      packer,
#                                      save_trajectories = TRUE)
#  
#  vcv <- diag(length(pars_list$start_pars)) * 0.00004
#  
#  sampler <- monty_sampler_random_walk(vcv)
# 
#  posterior <- likelihood + prior
# 
#   
#  samples <- monty_sample(posterior,
#                           sampler,
#                           200,
#                           initial = packer$pack(start_list),
#                           n_chains = 3)
#  
#  
# 
# qs_save(samples, outfile)

samples<-qs_read(outfile)

samples_df <- posterior::as_draws_df(samples)

summary<-posterior::summarise_draws(samples_df)

likelihood <- dust_likelihood_monty(filter_det, 
                                    packer,
                                    save_trajectories = TRUE)

posterior <- likelihood + prior

covmatrix<-cov(samples_df[,1:length(start_list)])
 
vcv<-2.38^2 / nrow(covmatrix) * covmatrix#

sampler <- monty::monty_sampler_adaptive(covmatrix)

#runner<-monty::monty_runner_callr(2)
runner <- monty::monty_runner_parallel(4)

samples <- monty::monty_sample(model=posterior, 
                               sampler=sampler, 
                               n_steps = nsamples,
                               initial = samples,
                               n_chains = nchains, 
                               runner = runner,
                               burnin = to_burn,
                               thinning_factor = each_n)


matplot(samples$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")



samples_df <- posterior::as_draws_df(samples)
bayesplot::mcmc_trace(samples_df)

summary<-posterior::summarise_draws(samples_df)



qs_save(samples, outfile)
qs_save(index,outfile_index)




