# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile         <-file.path(root, "output", "scenario_results.qs2")
outfile        <- file.path(root, "output", "scenario_plots.png")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(ggplot2)
library(dplyr)
library(data.table)



# Functions ---------------------------------------------------------------



get_scenario_quantiles<- function(dframe,varname){

  scen_names<-unique(dframe$scenario)  
  n=length(scen_names)
 
  datalist = vector("list", length = n)
  
  for (i in 1:n) {
    # ... make some data
    dat <- dframe[dframe$scenario==scen_names[i],] %>% group_by(days) %>% 
      do(data.frame(t(quantile(.[[varname]], probs = c(0.025,0.5,0.975)))))

    dat$scenario <- scen_names[i]  # maybe you want to keep track of which iteration produced it?
    
    datalist[[i]] <- dat # add it to your list
  }
  
  big_data = do.call(rbind, datalist)
  
}


# Load data ---------------------------------------------------------------


runs <- qs_read(infile)

my_names = c("days", 
             "scenario", 
             "cases_averted",
             "dose_per_cases_averted",
             "new_weekly_deaths",
             "new_yearly_hosp")

shortlist = lapply(runs, "[", , my_names)
shortDF <- rbindlist(shortlist)
shortDF$iteration<- rep(seq(1,500), each = 3651*4)

scen_tags<-c("Baseline",
             "Routine under 1yr",
             "Routine under 1 + campaign under 5yr",
             "Routine under 1 & 65yr")

shortDF$scenario<-factor(shortDF$scenario, levels = scen_tags)



# Prepare individual outputs ----------------------------------------------



# Cases averted

cases_averted<-get_scenario_quantiles(shortDF,"cases_averted") %>% 
  filter(scenario!="Baseline")%>% 
  filter(days==730 | days== 1825 | days==3650 )%>%
  mutate(year = case_when(days == 730 ~ 2, days == 1825 ~ 5 , days == 3650 ~ 10))

cases_averted$year<-factor(cases_averted$year, 
                                  levels=sort(unique(cases_averted$year)))

cases_averted$output<-"cases_averted"

# Dose per Case averted

dose_percase_averted<-get_scenario_quantiles(shortDF,"dose_per_cases_averted") %>% 
  filter(scenario!="Baseline") %>% 
  filter(days==730 | days== 1825 | days==3650 )%>%
  mutate(year = case_when(days == 730 ~ 2, days == 1825 ~ 5 , days == 3650 ~ 10))

dose_percase_averted$year<-factor(dose_percase_averted$year, 
                                  levels=sort(unique(dose_percase_averted$year)))

dose_percase_averted$output<-"dose_percase_averted"

# Deaths averted


df<- shortDF %>% 
  group_by(scenario,iteration) %>% 
  slice(seq(0, n(), by = 7)) %>% 
  mutate(cum_deaths=cumsum(new_weekly_deaths)) 

df <- df %>% 
  group_by(days,iteration) %>% 
  mutate(deaths_averted = cum_deaths[scenario == 'Baseline'] - cum_deaths )


deaths_averted<-get_scenario_quantiles(df,"deaths_averted") %>% 
  filter(scenario!="Baseline")%>% 
  group_by(scenario) %>% 
  mutate(week=seq(1,length(days))) %>% 
  filter(week==104 | week== 260 | week==520 )%>%
  mutate(year = case_when(week == 104 ~ 2, week == 260 ~ 5 , week == 520 ~ 10)) %>% 
  mutate(year = factor(year, levels=sort(unique(year))))

deaths_averted$output<-"deaths_averted"


# Hospital admissions


df<- shortDF %>% 
  group_by(scenario,iteration) %>% 
  slice(seq(0, n(), by = 365)) %>% 
  mutate(cum_admissions=cumsum(new_yearly_hosp)) 

df <- df %>% 
  group_by(days,iteration) %>% 
  mutate(admissions_averted = cum_admissions[scenario == 'Baseline'] - cum_admissions)


admissions_averted<-get_scenario_quantiles(df,"admissions_averted") %>% 
  filter(scenario!="Baseline")%>% 
  group_by(scenario) %>% 
  mutate(year=seq(1,length(days))) %>% 
  filter(year==2 | year== 5 | year== 10 ) %>% 
  mutate(year = factor(year, levels=sort(unique(year))))

admissions_averted$output<-"admissions_averted"


# Join dataframes

df<-rbind(cases_averted,dose_percase_averted,deaths_averted,admissions_averted)
df<-df %>% 
  mutate(output = case_when(
    output=="cases_averted" ~ "Norovirus AGE cases averted",
    output=="dose_percase_averted" ~ "Vaccine courses per case averted",
    output=="deaths_averted" ~ "Norovirus related deaths averted",
    output=="admissions_averted" ~ "Hospital admissions averted")) %>%
  mutate(output=factor(output,levels = unique(output)))

# Plot --------------------------------------------------------------------

options(scipen = 999) # turn off scientific notation

p<-ggplot(df, aes(x = year , y = X50.,color = scenario)) +
  geom_col(position = position_dodge(0.8), width = 0.7, fill = "white", alpha=0) +
  geom_errorbar(
    aes(ymin = X2.5., ymax = X97.5.),
    position = position_dodge(0.8), width = 0.2
  )+
  geom_point(aes(color = scenario), position = position_dodge(0.8)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "grey34"))+
  labs(x="Year",y=" ")+
  facet_wrap(~output,scales = "free") +
  theme_classic()+
  theme(
    legend.title = element_blank(),
    legend.position = "top",
  )



windows()
gridExtra::grid.arrange(p)

# save plot ---------------------------------------------------------------


ggsave(p, 
       filename = outfile ,
       device = "png")
