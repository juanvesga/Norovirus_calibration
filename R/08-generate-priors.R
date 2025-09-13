# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
outfile <- file.path(root, "output", "priors.qs2")
outfile2 <- file.path(root, "output", "limits.qs2")

outfile3 <- file.path(root, "output", "priors2.qs2")
outfile4 <- file.path(root, "output", "limits2.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(monty, include.only = c("monty_dsl"))

# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------

# mean
last_best=
  c(
    beta_1      =    0.13410824,
    beta_2      =    0.12708480,
    beta_3      =    0.20300807,
    beta_4      =    0.27545582,
    aduRR       =    0.58903565,
    maternalAB  =  367.16148967,
    imm_yr      =     16.15967444 ,
    imm_fac     =     1.70697752,
    repfac_0    =     64.26091774, 
    repfac_5    =     450.23810289,
    repfac_15   =     345.14095794,  
    repfac_65p  =    59.89184305,
  #  reported_var=      0.5,
    season_lag  =      30,
    season_amp  =      0.15
  )

prior <- monty_dsl({
  beta_1 ~     Uniform(0,1)#Gamma(shape = 2, scale = 0.15/2)
  beta_2 ~     Uniform(0,1)#Gamma(shape = 2, scale = 0.15/2)
  beta_3 ~     Uniform(0,1)#Gamma(shape = 3, scale = 0.2/3)
  beta_4 ~     Uniform(0,1)#Gamma(shape = 3, scale = 0.2/3)
  aduRR  ~     Uniform(0,1)
  maternalAB ~ Gamma(shape = 3, scale = 150/3)
  imm_yr~      Uniform(0,100)#Gamma(shape = 4, scale = 20/4)
  imm_fac~     Uniform(0,5)
  repfac_0 ~   Uniform(0,900)#Gamma(shape = 6, scale = 128/6)
  repfac_5 ~   Uniform(0,1500)#Gamma(shape = 12, scale = 523/12)
  repfac_15~   Uniform(0,900)#Gamma(shape = 12, scale = 500/12)
  repfac_65p ~ Uniform(0,900)#Gamma(shape = 6, scale = 98/6)
#  reported_var~ Uniform(0,1)
  season_lag  ~ Uniform(-10+20,15+20)
  season_amp  ~ Uniform(0,1)
})

prior2 <- monty_dsl({
  beta_1 ~     Uniform(0,1)#Gamma(shape = 2, scale = 0.15/2)
  beta_2 ~     Uniform(0,1)#Gamma(shape = 2, scale = 0.15/2)
  beta_3 ~     Uniform(0,1)#Gamma(shape = 3, scale = 0.2/3)
  beta_4 ~     Uniform(0,1)#Gamma(shape = 3, scale = 0.2/3)
  aduRR  ~     Uniform(0,1)
  maternalAB ~ Gamma(shape = 3, scale = 150/3)
  imm_yr~      Uniform(0,100)#Gamma(shape = 4, scale = 20/4)
  repfac_0 ~   Uniform(0,900)#Gamma(shape = 6, scale = 128/6)
  repfac_5 ~   Uniform(0,1500)#Gamma(shape = 12, scale = 523/12)
  repfac_15~   Uniform(0,900)#Gamma(shape = 12, scale = 500/12)
  repfac_65p ~ Uniform(0,900)#Gamma(shape = 6, scale = 98/6)
#  reported_var~ Uniform(0,1)
  season_lag  ~ Uniform(-10+20,15+20)
  season_amp  ~ Uniform(0,1)
})



# mean
limits<-list(
lower=
  log(c(
    beta_1      =    0,
    beta_2      =    0,
    beta_3      =    0,
    beta_4      =    0,
    aduRR       =    0,
    maternalAB  =    7,
    imm_yr      =     0 ,
    imm_fac     =     0,
    repfac_0    =     0, 
    repfac_5    =     0,
    repfac_15   =     0,  
    repfac_65p  =     0,
 #   reported_var=     0,
    season_lag  =     -10+20,
    season_amp  =     0
    
  )),
upper=
  log(c(
    beta_1      =     1,
    beta_2      =     1,
    beta_3      =     1,
    beta_4      =     1,
    aduRR       =     1,
    maternalAB  =     700,
    imm_yr      =     100,
    imm_fac     =     10,
    repfac_0    =     900, 
    repfac_5    =     1500,
    repfac_15   =     900,  
    repfac_65p  =     900,
 #   reported_var=     0.5,
    season_lag  =     15+20,
    season_amp  =     1
  ))
)


limits2<-list(
  lower=
    log(c(
      beta_1      =    0,
      beta_2      =    0,
      beta_3      =    0,
      beta_4      =    0,
      aduRR       =    0,
      maternalAB  =    7,
      imm_yr      =     0 ,
      repfac_0    =     0, 
      repfac_5    =     0,
      repfac_15   =     0,  
      repfac_65p  =     0,
   #   reported_var=     0,
      season_lag  =     -10+20,
      season_amp  =     0
      
    )),
  upper=
    log(c(
      beta_1      =     1,
      beta_2      =     1,
      beta_3      =     1,
      beta_4      =     1,
      aduRR       =     1,
      maternalAB  =     700,
      imm_yr      =     100,
      repfac_0    =     900, 
      repfac_5    =     1500,
      repfac_15   =     900,  
      repfac_65p  =     900,
 #     reported_var=     0.5,
      season_lag  =     15+20,
      season_amp  =     1
    ))
)


# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(prior, outfile)
qs_save(limits, outfile2)
qs_save(prior2, outfile3)
qs_save(limits2, outfile4)
