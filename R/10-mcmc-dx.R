# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")

outfile1        <- file.path(root, "output", "mcmc_fits.qs2")
outfile2        <- file.path(root, "output", "mcmc_inits.qs2")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))
source(file.path(root, "R", "collect_function.R"))
# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(abind, include.only = c("abind"))

library(ggplot2)
library(dplyr)
library("lattice")  ## for the 'xyplot' command
library(fitR)
# -------------------------------------------------------------------------
#  -------------------------------------------------------------------
# -------------------------------------------------------------------------

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

burn_in    <- 0
thin_in    <- 30


infile_chains <- file.path(root, "output", paste0("mcmc_chains",model,".qs2"))


trace     <- qs_read(infile_chains) 
d<- ncol(trace[[1]])
trace[[1]][,1:(d-1)] <-exp(trace[[1]][,1:(d-1)])
trace[[2]][,1:(d-1)] <-exp(trace[[2]][,1:(d-1)])
trace[[3]][,1:(d-1)] <-exp(trace[[3]][,1:(d-1)])
trace[[4]][,1:(d-1)] <-exp(trace[[4]][,1:(d-1)])


# burn and thin
traceBurn <- burnAndThin(trace, burn = burn_in)
chains <- burnAndThin(traceBurn, thin = thin_in)



plot1<-xyplot(chains,layout = c(4, 4))
gridExtra::grid.arrange(plot1)


# removing the burn-in increases the ESS
print(coda::effectiveSize(chains))


##           R_0         D_lat         D_inf         alpha         D_imm 
##      2495.223      3024.481      3132.444      3201.725      3188.532 
##           rho      logPrior logLikelihood    logDensity 
##      3350.921         0.000      2681.293      2681.293

# autocorrelation

plot3<-coda::acfplot(chains, lag.max = 60)
gridExtra::grid.arrange(plot3)



# Note that plotPosteriorDensity can take a list of mcmc.list It will plot the
# different mcmc.list by combining their elements Let's plot the combined
# unthinned trace vs the combined thinned trace.
plot4<-plotPosteriorDensity(list(full =traceBurn , thinned_burn = chains))
gridExtra::grid.arrange(plot4)

plot5<-levelplot(chains[[1]], col.regions = heat.colors(100))
gridExtra::grid.arrange(plot5)


post<-rbind(  chains[[1]][,"logDensity"],
              chains[[2]][,"logDensity"],
              chains[[3]][,"logDensity"],
              chains[[4]][,"logDensity"])

matplot(t(post),type = "l")


