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
infile_chains  <- file.path(root, "output", "mcmc_chains.qs2")
infile_inits    <- file.path(root, "output", "mcmc_inits.qs2")
#infile_index   <- file.path(root, "output", "index.qs")
infile_odin1chain   <- file.path(root, "output", "chains_imm2par.RData")

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


# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------


# Source Odin model
source(infile_model2)
parameters <- qs_read(infile0)
pars_list  <- qs_read(infile_input)
chains     <- qs_read(infile_chains) 
inits      <- qs_read(infile_inits) 
load(infile_odin1chain) 
chains_odin <-processed_chains$pars
rm(processed_chains)


#combine chains
longchain<-do.call(rbind, chains)



#get last sample
#params<- unlist(exp(longchain[nrow(longchain),1:(ncol(longchain)-1)]))

params<- unlist(exp(chains_odin[nrow(chains_odin),]))

# Set staring parameters
pars<-update_parameters(params,pars_list$pars_list)  

#dust object
#sys  <- dust_system_create(noro_model(), pars, n_particles = 1)

# Determinsitic object
sys  <- dust_system_create(noro_model(), pars, deterministic = TRUE, n_particles = 1)

# Test runs ---------------------------------------------------------------


# # Get state index
# index<- dust_unpack_index(sys)
# 
# #Get default state initial
state<-dust_system_state(sys)
# 
named_state<-dust_unpack_state(sys,state)
# 
# # Set inital model state based on calibrated model
# 
# # Set default initial state of the model 
dust_system_set_state_initial(sys)
# 
# # Set parameters
dust_system_update_pars(sys, pars=pars)
# 
# 
# run simulation
endsim<-24109
tt<-seq(0,endsim)
y0<- dust_system_simulate(sys, tt)
# 
# # unpack results
y0 <- dust_unpack_state(sys, y0)
# 
# 
# # Get model state at the latest point for future simulations
# init_state0<-dust_system_state(sys)
# 
# names(init_state0)<-paste0(names(y0))
# 
# 
cumcases<-(cumsum(y0$new_cases))
# 
plot(tt,(y0$new_cases),
        type="l",
        xlim = c(2,30000),
        ylim = c(0, 1e5))
# 
# 
plot(tt[tt %% 7 == 0], (y0$new_cases_week[tt %% 7 == 0]), type = "l",
     ylab = "Weekly incidence",
     xlab = "Time",
     xlim = c(15e3,25e3),
     ylim=c(0,4e5))
# 
plot(tt[tt %% 7 == 0], y0$new_cases_week_gii4[tt %% 7 == 0], type = "l",
     ylab = "Weekly incidence",
     xlab = "Time",
     xlim = c(2000,10000),
     ylim=c(0,1e5))
# 
# 
rate<-100/1e5
plot(tt[tt %% 7 == 0], rate*(y0$new_cases_week[tt %% 7 == 0]), type = "l",
     ylab = "Weekly admissions",
     xlab = "Time",
     xlim = c(15e3,25e3),
     ylim=c(0,250)
)


plot(tt[tt %% 7 == 0], colSums(y0$reported_wk[,tt %% 7 == 0]), type = "l",
     ylab = "Weekly reported",
     xlab = "Time",
     xlim = c(15e3,25e3),
     ylim=c(0,1000)
)






### Interventions

init_state0<-unlist(inits[nrow(inits),])

# Primer run for 5 yrs
t_zero           <- 24110
t_span           <-365*5
pars0<-list() # no updates

runs_primer<-update_interventions(sys,t_zero,init_state0,t_span,pars0,"Baseline")


# Get model state at the latest point for future simulations
init_state0<-dust_system_state(sys)

#names(init_state0)<-paste0(names(y0))

# 
# 
# 
# 
# 
# 
# 
# # Get state index
# index<- dust_unpack_index(sys)
# 
# #Get default state initial
# state<-dust_system_state(sys)
# 
# named_state<-dust_unpack_state(sys,state)
# 
# # Set inital model state based on calibrated model
# 
# # Set default initial state of the model 
# dust_system_set_state_initial(sys)
# 
# # Set parameters
# dust_system_update_pars(sys, pars=pars)
# 
# 
# # run simulation
# endsim<-24109
# tt<-seq(0,endsim)
# y0<- dust_system_simulate(sys, tt)
# 
# # unpack results
# y0 <- dust_unpack_state(sys, y0)
# 
# 
# # Get model state at the latest point for future simulations
# init_state0<-dust_system_state(sys)
# 
# names(init_state0)<-paste0(names(y0))
# 
# 
# cumcases<-(cumsum(y0$new_cases))
# 
# plot(tt,(y0$new_cases),
#         type="l",
#         xlim = c(2,30000),
#         ylim = c(0, 1e5))
# 
# 
# plot(tt[tt %% 7 == 0], (y0$new_cases_week[tt %% 7 == 0]), type = "l",
#      ylab = "Weekly incidence", 
#      xlab = "Time",
#      xlim = c(15e3,25e3),
#      ylim=c(0,4e5))
# 
# plot(tt[tt %% 7 == 0], y0$new_cases_week_gii4[tt %% 7 == 0], type = "l",
#      ylab = "Weekly incidence", 
#      xlab = "Time",
#      #xlim = c(2000,10000),
#      ylim=c(0,1e5))
# 
# 
# rate<-100/1e5
# plot(tt[tt %% 7 == 0], rate*(y0$new_cases_week[tt %% 7 == 0]), type = "l",
#      ylab = "Weekly admissions", 
#      xlab = "Time",
#      xlim = c(15e3,25e3),
#      ylim=c(0,250)
# )

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


vacc_eff         <- 0.75
vacc_coverage    <- 0.9
campaign_coverage<- 0.7
t_zero           <-1
t_span           <-365*10
t_campaign_stop  <-t_zero+180
vacc_dis         <- 1
vacc_trans       <- 0



# Baseline

pars0<-pars

runs0<-update_interventions(sys,t_zero,init_state0,t_span,pars0,"Baseline")
runs0$cases_averted<-runs0$cum_new_daily_cases-runs0$cum_new_daily_cases
runs0$dose_per_cases_averted<-runs0$cum_vacc_doses/(runs0$cum_new_daily_cases-runs0$cum_new_daily_cases)


deaths<- as.data.frame(runs0$new_weekly_deaths) %>% 
  slice(seq(0, n(), by = 7))

plot(deaths$`runs0$new_weekly_deaths`, type = "l", col="black",
     ylab = "Weekly Noro deaths", 
     xlab = "Week"
     #xlim = c(2000,10000),
     #ylim=c(0,1e5)
     )


adm<-data.frame(runs0$new_yearly_hosp)

admissions<- as.data.frame(runs0$new_yearly_hosp) %>% 
  slice(seq(0, n(), by = 365))

plot(admissions$`runs0$new_yearly_hosp`, type = "l", col="black",
     ylab = "Yearly Norovirus Admissions ", 
     xlab = "year"
)


reported<-as.data.frame(runs0$new_weekly_reported) %>% 
  slice(seq(0, n(), by = 7))

plot(reported$`runs0$new_weekly_reported`, type = "l", col="black",
     ylab = "Weekly reported ", 
     xlab = "Week"
)


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








# Join together all results 
#############  by strain
scen_tags<-c("Baseline",
             "Routine under 1yr",
             "Routine under 1 & 65yr",
             "Routine under 1 + campaign under 5yr")


df<-rbind(runs0,runs1,runs2,runs3)
df$scenario<- factor(df$scenario, levels=scen_tags)

df2<- data.frame(
  d0=runs0$new_weekly_deaths,
  d1=runs1$new_weekly_deaths,
  d2=runs2$new_weekly_deaths,
  d3=runs3$new_weekly_deaths)

dfdeaths<- df2  %>%
  slice(seq(0, n(), by = 7))

dfdeaths<-cumsum(dfdeaths)

dfdeaths<-dfdeaths[nrow(dfdeaths),]

averted_deaths<-dfdeaths$d0-dfdeaths[2:4]


df2<- data.frame(
  d1=runs1$cases_averted,
  d2=runs2$cases_averted,
  d3=runs3$cases_averted)

averted_cases<-df2[nrow(df2),]






plot1<-ggplot(df, aes(x=days,y=cum_new_daily_cases,colour = scenario))+
  geom_line()+
  labs(title = "Vaccination effect on Norovirus AGE cases", y="Cumulative daily Incidence")+
  #ylim(0,75)+
  theme_minimal()




plot2<-ggplot(df, aes(x=days,y=new_daily_cases,colour = scenario))+
  geom_line()+
  labs(title = "Vaccination effect on Norovirus AGE cases", y="Daily Incidence")+
  #ylim(0,75)+
  theme_minimal()

windows()
gridExtra::grid.arrange(plot1,plot2)


# dose per case averted
df2<-df %>% filter(days==t_span) %>% filter(scenario!="Baseline") %>% 
  select(c("dose_per_cases_averted","scenario"))



plot3<-ggplot(df2, aes(x=scenario,y=dose_per_cases_averted, fill=scenario))+
  geom_bar(stat = "identity")+
  labs(x="Scenario",y="Dose per case averted",title = "Dose per cases averted after 5 years")+
  theme(legend.position="none")


# cases averted 
tmp<-data.frame(
  scenario= factor(scen_tags[2:4],levels = scen_tags[2:4]),
  cases_averted=as.numeric(averted_cases[1,]))
    
    
plot4<-ggplot(tmp, aes(x=scenario,y=cases_averted, fill=scenario))+
  geom_bar(stat = "identity")+
  labs(x="Scenario",y="Age cases averted",title = " Cases averted after 10 years")+
  theme(legend.position="none")



# Early Deaths prevented 
tmp<-data.frame(
  scenario= factor(scen_tags[2:4],levels = scen_tags[2:4]),
  deaths_averted=as.numeric(averted_deaths[1,]))


plot5<-ggplot(tmp, aes(x=scenario,y=deaths_averted, fill=scenario))+
  geom_bar(stat = "identity")+
  labs(x="Scenario",y="Age related deaths averted",title = " Deaths averted after 10 years")+
  theme(legend.position="none")






windows()
gridExtra::grid.arrange(plot3,plot4,plot5, nrow=2)






get_scen_cases<-function(scen,scen_text){
  
  scen<-scen %>% filter(
    scenario==scen_text
  )
  
  
  df<-data.frame(c(
    scen$new_weekly_cases_gi3, # GI3
    scen$new_weekly_cases_gi,  # Other Gi
    scen$new_weekly_cases_gii4,
    scen$new_weekly_cases_gii
  ))
  
  df$strain<-factor(c(
    rep("GI3",length(scen$new_weekly_cases_gi3)),
    rep("GI",length(scen$new_weekly_cases_gi3)),
    rep("GII4",length(scen$new_weekly_cases_gi3)),
    rep("GII",length(scen$new_weekly_cases_gi3))
  ), levels =c("GI3","GI","GII4","GII"))
  

  
  df$date<-c(
    rep(scen$days,4))
  
  names(df)<-paste(c("cases","strain","date"))
  
  return(df)
}

df1<-get_scen_cases(df,scen_tags[1])
df2<-get_scen_cases(df,scen_tags[2])
df3<-get_scen_cases(df,scen_tags[3])
df4<-get_scen_cases(df,scen_tags[4])


get_scen_plot<-function(df,tag){  
  

  data<- df %>% group_by(strain) %>%
    slice(seq(1, n(), by = 7))
  
  data$week<-rep(seq(1,length(unique(data$date))),4)
  
  gp<-ggplot(data, aes(fill=strain, y=cases, x=week,width=1)) + 
    geom_bar(position="stack",stat="identity")+
    scale_fill_manual(values=c('skyblue3','firebrick2','yellow2','grey28'))+
    # scale_x_date(date_breaks = "4 months", date_labels = "%b",
    #              limits = as.Date(c('2025-12-01','2029-08-31'))) +
    labs(y="AGE cases",x="Week",title=tag)+
    ylim(c(0,5e4))+
    xlim(c(0,52*4))+
    theme_minimal() +  
    theme(
      #legend.position = "none",
      title.text = element_text(size = 10),
      panel.background = element_blank(),
      axis.text = element_text(colour = "black", size = 10),
      axis.title.y = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9), legend.key = element_blank(),
      #axis.text.x = element_text(angle = 60, hjust = 1, size=8),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0.1, 0.95),
      panel.border = element_rect(colour = "black", fill=NA)
    )
  
  return(gp)
  
}

p0<-get_scen_plot(df1,scen_tags[1])
p1<-get_scen_plot(df2,scen_tags[2])
p2<-get_scen_plot(df3,scen_tags[3])
p3<-get_scen_plot(df4,scen_tags[4])

windows()
gridExtra::grid.arrange( p0,p1,p2,p3,
                         layout_matrix = rbind(c(1, 2),
                                               c(3, 4)))

