# Check packages ----------------------------------------------------------


required <- c(
  "here",
  "qs2",
  "readxl",
  "data.table",
  "gridExtra",
  "odin2",
  "dust2",
  "monty",  
  "tools",
  "matrixStats",
  "lubridate",
  "ISOweek",
  "tidyr",
  "dplyr",
  "ggplot2")

installed <- rownames(installed.packages())
not_installed <- required[!required %in% installed]

if (length(not_installed)) {
  msg <- "Some required packages are missing. Please install the following:\n"
  pkgs <- paste("\t", not_installed, collapse = "\n")
  msg <- sprintf("%s%s", msg, pkgs)
  stop(msg)
} else {
  # Set the root directory --------------------------------------------------
  root <- here::here()
  
  # List the fils we need to run in order
  files <- file.path(
    root, "R",
    c(
      "01-parameters.R",
      "02-data-calibration.R",
      "03-data-polymod.R",
      "04-data-school-uk.R",
      "05-data-comix.R",
      "06-data-covid.R",
      "07-generate-input-parameters.R",
      "08-generate-priors.R",
      "09-mcmc.R",
      "10-mcmc-dx.R"

    )
  )
  
  # Loop over the files
  for (f in files) {
    local(source(f, local = TRUE))
  }
  
  message("Analysis complete. Results are in the 'results/' directory.")
}



# 
# library(odin2)
# library(dust2)
# library(monty)
# 
# root<-here::here()
# path<-here::here("src","model_2pars.R")
# 
# #odin2::odin_migrate(path = path, dest = path)
# # Call model object
# source(path)
# 
# # Call data 
# source(here::here("src","utility_functions.R"))
# source(here::here("scripts","load_data.R"))
# source(here::here("scripts","setup.model.R"))
# 
# 
# sys <- dust_system_create(noro_model(), pars = list(), dt = 1)
# dust_system_set_state_initial(sys)
# 
# 
# 
# incidence<- c(2, 4, 1, 3, 2, 5)
# cases    <- c(12,23,25,36,30,57)
# 
# dpois(cases/600, incidence/600, log = TRUE)

