# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
outfile <- file.path(root, "output", "parameters.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
modify_attach(qs2, include.only = "qs_save")


# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------

tstep <- "day"

if(tstep=="day"){
  t_step=1
}else
{
  t_step=7
}

sim_startdate= epi_week("1960-01-01")
sim_enddate  = epi_week("2025-12-31")

parameters<-list(
  tstep=tstep,
  t_step=t_step,
  sim_startdate=sim_startdate,
  sim_enddate  = sim_enddate, 
  time_vec     = seq(sim_startdate,sim_enddate,tstep),
  weeks_vec    = seq(sim_startdate,sim_enddate,"week"),
  years_vec    = seq(sim_startdate,sim_enddate,"year"),
  mort_rates =c(   
    0.0039640, #1
    0.0002937, #2
    0.0001537, #3
    0.0001155, #4
    0.0000924, #5
    0.0000924, #6
    0.0000874, #7
    0.0000891, #7-15
    0.0002980, #16-25 
    0.0005477, #26-35
    0.0011747, #36-45
    0.0026180, #46-55
    0.0064792, #56-65
    0.0163398+(1/(90-75))), #66-75 ##<--- adding the aging from 75 to 100 as a limit 
    ages   =  c(1,2,3,4,5,6,7,15,25,35,45,55,65,75), # Upper end of age bands
  vaccination_coverage = seq(1,14)*0,
  cfr        =c(
    0.006658405,
    0.006930603,
    0.007165995,
    0.00736015,
    0.007712279,
    0.00797422,
    0.008300209,
    0.009486607,
    0.013072195,
    0.018749172,
    0.026713006,
    0.038570177,
    0.053864667,
    0.091299783),
  crossp_GI =0.05,
  crossp_GII =0.05
)



# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(parameters, outfile)
