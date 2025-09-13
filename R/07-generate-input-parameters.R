
# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile1  <- file.path(root,"output", "parameters.qs2")
infile2  <- file.path(root,"output", "polymod.qs2")
infile3  <- file.path(root,"output", "comix.qs2")
infile4  <- file.path(root,"output", "school_uk.qs2")
infile5  <- file.path(root,"output", "covid_sche.qs2")

outfile <- file.path(root, "output", "params_list.qs2")
#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(monty, include.only = c("monty_packer"))

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------

parameters <- qs_read(infile1)
cmat       <- qs_read(infile2)
comix      <- qs_read(infile3)
school_uk  <- qs_read(infile4)
covid_sche <- qs_read(infile5)


fixed_pars<-list(
  pop=parameters$pop,
  mu =parameters$mort_rates/(365),
  cfr=parameters$cfr,
  aging_vec = c(1/head(diff(c(0,parameters$ages)),-1),0)/(365), 
  N_age =length(parameters$ages),
  m = cmat$transmission, 
  m_holi=cmat$transmission_holi, 
  cmx_1 =comix$cmx_1,
  cmx_2 =comix$cmx_2,
  cmx_3 =comix$cmx_3,
  cmx_4 =comix$cmx_4,
  cmx_5 =comix$cmx_5,
  cmx_6 =comix$cmx_6,
  cmx_7 =comix$cmx_7,
  cmx_8 =comix$cmx_8,
  cmx_9 =comix$cmx_9,
  school_time = seq(0,length(school_uk)-1),
  school_value= school_uk,
  comix_time  = seq(0,length(covid_sche)-1),
  comix_value = covid_sche,
  vaccination_coverage=parameters$vaccination_coverage,
  vacc_switch_on = seq(1,length(parameters$ages))*0,
  campaign_switch= seq(1,length(parameters$ages))*0,
  vacc_trans = 0,
  vacc_dis = 0,
  crossp_GI =parameters$crossp_GI,
  crossp_GII =parameters$crossp_GII
)

pars<-c("beta_1",
        "beta_2",
        "beta_3",
        "beta_4",
        "aduRR",
        "maternalAB",
        "imm_yr", 
        "imm_fac",
        "repfac_0",
        "repfac_5",
        "repfac_15",
        "repfac_65p")

packer<-monty_packer(pars, fixed=fixed_pars)



start_pars<-list(
  beta_1=0.13953720,
  beta_2=0.12680114,
  beta_3=0.20410557,
  beta_4=0.27306245,
  aduRR=0.69372532  ,
  maternalAB=520.61197156,
  imm_yr=17.18686356,
  imm_fac=1.68360017 ,
  repfac_0=105.07189709,
  repfac_5=954.17803464 ,
  repfac_15=581.94243972 ,
  repfac_65p=95.41216406  ,
  reported_var=0.4
)




# start_pars<-list(
#   beta_1      =    0.20752229,   
#   beta_2      =    0.26123320,   
#   beta_3      =    0.21955518,   
#   beta_4      =    0.18338303,   
#   aduRR       =    0.66758771, 
#   maternalAB  =  154.10159618,   
#   imm_yr      =     5.39353643,   
#   imm_fac     =     2.33102585,  
#   repfac_0    =     54.82729580, 
#   repfac_5    =     379.84666082, 
#   repfac_15   =     389.28583888 ,   
#   repfac_65p  =    68.25388766,   
#   crossp_GI   =        0.15067537,   
#   crossp_GII  =      0.05961897)



# Create filter
input_object<-list(
start_pars=start_pars,
pars_list=c(start_pars,fixed_pars),
packer   = packer,
free_pars_names=pars,
fixed_pars=fixed_pars)

# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(input_object, outfile)



