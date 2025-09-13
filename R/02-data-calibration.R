# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile0  <- file.path(root,"output", "parameters.qs2")
infile1  <- file.path(root,"data", "raw", "IID2_incidence_OBrien.csv")
infile2  <- file.path(root,"data", "raw", "sgss_all_age.csv")
infile3  <- file.path(root,"data", "raw", "sgss_strain_individual.csv")
infile4  <- file.path(root,"data", "raw", "serology_prev.csv")
infile5  <- file.path(root,"data", "raw", "admissions_sandman.csv")

outfile1 <- file.path(root, "output", "data_long.qs2")
outfile2 <- file.path(root, "output", "data_short.qs2")
outfile3 <- file.path(root, "output", "data_for_plots.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(lubridate, include.only = c("ceiling_date","isoweek","ymd","weeks","years"))
modify_attach(ISOweek, include.only = c("ISOweek2date"))
modify_attach(imputeTS,include.only = c("na_interpolation"))

require(dplyr)

# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------


parameters <- qs_read(infile0)


tstep <- parameters$tstep
t_step<- parameters$t_step

time_vec<-parameters$time_vec


# IID2 incidence data-------------------------------------------------------

# Source: IID2 survey (published )
# Type: Community incidence per 1000 person-years
# Strata: Age structured (five age groups)
# Time span: one data point in 2008

idd2_startdate<-epi_week("2008-01-01")
idd2_enddate<-epi_week("2008-12-28")
idd2_fup<- seq(idd2_startdate,idd2_enddate,tstep)
index_idd2<-which(time_vec %in% idd2_startdate)

st<-which(time_vec==idd2_startdate)
ed<-which(time_vec==idd2_enddate)
span<-seq(st,ed,1) # time span of idd2 in the sim
iid2_time<-span[ which(span%%t_step ==0) ]
iid2_time<-iid2_time[floor(length(iid2_time)/2)] # for day average take

raw<- read.csv(infile1, header=TRUE)
raw$ci<-as.numeric(raw$CI_upper-raw$CI_lower)
raw$CI_lower<-as.numeric(raw$CI_lower)
raw$CI_upper<-as.numeric(raw$CI_upper)


data_iid2.c4<-raw %>%
  filter (age_cat!="5-") 

data_iid2.c4$age_cat<- factor(data_iid2.c4$age_cat, 
                              levels = c("0 to 1","1 to 5","5 to 15","15 to 64","65+"))
# Keep for plotting  
data_iid2.c4<-data_iid2.c4[order(data_iid2.c4$age_cat),]

#data object
df_iid2<-data.frame(time=iid2_time,
                    cases_a1=data_iid2.c4$per1000personyears[1],
                    cases_a1=data_iid2.c4$per1000personyears[2],
                    cases_a1=data_iid2.c4$per1000personyears[3],
                    cases_a1=data_iid2.c4$per1000personyears[4],
                    cases_a1=data_iid2.c4$per1000personyears[5])


# SGSS Weekly series all-------------------------------------------------------

# Source: UKHSA surveillance of clinical cases reported, and molecular character
# Type: Time-series of reported cases 
# Strata: aggregated in strain structrure and age  
# Time span: Weekly "2014-01-05" to "2023-07-02" but used until jan 2020

raw<- read.csv(infile2, header=TRUE)

startd<-paste(raw$Year[1],"-W0",raw$Reporting_week[1],"-1",sep="")
startd<-ISOweek2date(startd)
dates2<-epi_week(seq(startd, by = "week", length.out = nrow(raw)))
sgss_days<-which(time_vec%in%dates2)

df<-t(as.matrix(raw[,c("a0.4", "a5.15", "a15.64", "a65")]))
rownames(df) <- c("a0_4","a15_64","a5_14","a65p")


#### Total cases series
total_cases<-data.frame(cases=colSums(df, na.rm = T))
total_cases$day<-sgss_days
total_cases$date<-dates2
id<-which(total_cases$date<epi_week("2020-01-01")) # fit to before covid

# Keep for plotting
total_cases<-total_cases[id,]
total_cases$week<-as.factor(isoweek(total_cases$date))
# data object
df_sgss_all<-data.frame(
  time=total_cases$day,
  reported=total_cases$cases
)


# Get average seasons
  require(ISOweek)
  
  fiveyr_av<-total_cases %>% 
    group_by(week) %>% 
    summarize(average=round(mean(cases)))
  
  
  
  fiveyr_av$yweek<-paste("2016",fiveyr_av$week,"3",sep=":")
  
  fiveyr_av$date<-epi_week(as.Date(fiveyr_av$yweek, format = "%Y:%W:%u")) 
  fiveyr_av<-fiveyr_av[-c(nrow(fiveyr_av)-1,nrow(fiveyr_av) ),]
  fiveyr_av$day<-which(time_vec%in%fiveyr_av$date)
  
  half1<-fiveyr_av[c(1:30),]
  half2<-fiveyr_av[c(31:51),]
  
  fiveyr_av2<-rbind(half2,half1)
  fiveyr_av2$date<-seq(from=half2$date[1],by="week", length.out=51)
  fiveyr_av2$day<-which(time_vec%in%fiveyr_av2$date)
  
  ## Comment this out to get the full time series 
  total_cases_avg<-data.frame(
    day=fiveyr_av2$day, 
    cases=fiveyr_av2$average,
    date=fiveyr_av2$date)

  # data object
  df_sgss_avg<-data.frame(
    time=total_cases_avg$day,
    reported_avg=total_cases_avg$cases
  )


 

# SGSS proportion reported by age-----------------------------------------------

# Source: UKHSA surveillance of clinical cases reported, and molecular character
# Type: aggregated reported data by age
# Strata: By age  
# Time span: Weekly "2014-01-05" to "2023-07-02" but used until jan 2020

# Keep for plotting
agg_age<-(rowSums(df, na.rm = T)/sum(rowSums(df, na.rm = T)))%*%t(c(0.8, 1,1.2))


# weight population
a0_size<-56308/12
a5_size<-56308/12
a15_size<-56308/12
a65_size<-56308/12


# data object
df_sgss_age<-data.frame(
  time=df_sgss_all$time[nrow(df_sgss_all)],
  a0_event  =round(agg_age[1,2]*a0_size),
  a5_event  =round(agg_age[2,2]*a5_size),
  a15_event =round(agg_age[3,2]*a15_size),
  a65_event =round(agg_age[4,2]*a65_size))


# SGSS series by Strain-----------------------------------------------

# Source: UKHSA surveillance of clinical cases reported, and molecular character
# Type: Time-series of reported cases by strain
# Strata: By model strain structrure   
# Time span: Weekly 2018  to 2023 but used until jan 2020


raw<- read.csv(infile3, header=TRUE)

id1<-which(raw$Characterisation.results=="GI indeterminate")
raw<-raw[-id1,]
id2<-which(raw$Characterisation.results=="GII indeterminate")
raw<-raw[-id2,]
id3<-which(raw$Characterisation.results=="GII.4 and GI.3")
raw<-raw[-id3,]
id4<-which(raw$Characterisation.results=="GII.4 and other GI")
raw<-raw[-id4,]
id5<-which(raw$Characterisation.results=="Mixed indeterminate")
raw<-raw[-id5,]
id6<-which(raw$Characterisation.results=="Other GII and other GI")
raw<-raw[-id6,]
table(raw$Characterisation.results)


raw<-raw %>%
  group_by(Year,ISOweek,Characterisation.results) %>%
  summarize(n_observations = n())


raw_wide <- tidyr::spread(raw, Characterisation.results, n_observations)
head(raw_wide)
dates_wk<-ymd("2018-01-01") + weeks(raw_wide$ISOweek - 1) + 
  years(raw_wide$Year - raw_wide$Year[1])


dates<-seq(epi_week(dates_wk[1]),epi_week(dates_wk[length(dates_wk)]),tstep)


data<-t(as.matrix(raw_wide[,3:6]))
rownames(data) <- c("GI.3","GII.4","O-GI","O-GII")
id<-which(is.na(data[1,]))
data[1,id]<-0
id<-which(is.na(data[2,]))
data[2,id]<-0
id<-which(is.na(data[3,]))
data[3,id]<-0
id<-which(is.na(data[4,]))
data[4,id]<-0

# plot(dates_wk,data[2,],col="blue",type="o")
# lines(dates_wk,data[3,],col="red",type="l")
# lines(dates_wk,data[2,], col="gold",type="l")
# lines(dates_wk,data[4,],col="grey23",type="l")
# 
# rowSums(data)



df<-data.frame(x=dates_wk,
               week=isoweek(dates_wk),
               gi3=data[1,],
               gi=data[3,],
               gii4=data[2,],
               gii=data[4,])

# Stop data at start of covid-19

remove<-which(dates_wk=="2020-03-26")

df<-df[1:remove,]

df$day<-which(time_vec%in%df$x)

# Keep for plotting
total_cases_str<-df

# Data object
df_sgss_strain<-data.frame(
  time=total_cases_str$day,
  reported_gi3 =total_cases_str$gi3,
  reported_gi  =total_cases_str$gi,
  reported_gii4=total_cases_str$gii4,
  reported_gii =total_cases_str$gii
)





# GII4 seroprevalebce in children-----------------------------------------------

# Source: Cross sectional of serology among children in England (Lindesmith et al)
# Type: Prevalence data from community survey
# Strata: By age 0-7 
# Time span: "2011-07-01"


sero1_n<-103
sero2_n<-107
sero3_n<-121
sero4_n<-124
sero5_n<-122
sero6_n<-109

raw<- read.csv(infile4, header=TRUE)

names(raw)<-paste(c("age","mean","low","up"))

# Keep for plotting
dfsero<-raw

date_sero<-epi_week("2011-07-01")
day_sero<-which(time_vec%in%date_sero)

# data object
df_gii4_kids<-data.frame(
  time=day_sero,
  sero1=round(raw$mean[1]*sero1_n),
  sero2=round(raw$mean[2]*sero2_n),
  sero3=round(raw$mean[3]*sero3_n),
  sero4=round(raw$mean[4]*sero4_n),
  sero5=round(raw$mean[5]*sero5_n),
  sero6=round(raw$mean[6]*sero6_n))


# Hospital Admissions -------------------------------------------------------

# Source: Sandman et al for primary admissions and age distributions
# Type: Time-series of estimated admissions 
# Strata: added age structure   
# Time span: Weekly "2009-01-20" to "2016-06-01" 

raw<- read.csv(infile5, header=TRUE)

raw$eweekdate<- epi_week(raw$date) # format to epiweek

raw$diff<- raw$eweekdate-as.Date(raw$date) # find timediff for duplicates

# remove duplicates with the largest difference   
df0<- raw %>% 
  group_by(eweekdate) %>% 
  slice_min(diff, n=1) %>% 
  ungroup %>% 
  select(cases,eweekdate) %>% 
  rename(date = eweekdate)

# create master vector of dates with complete sequence
fullwk<-seq(ymd('2009-07-09'),ymd('2016-07-09'), by = '1 week')

week<-data.frame(date=epi_week(fullwk),week=lubridate::epiweek(fullwk))


# Merge data and complete dates
df<-left_join(week, df0, by=c("date")) %>% 
  select(c("date","week","cases"))

df$year<-lubridate::year(df$date)

# Interpolate on missing gaps 
df1  <- df %>%
  group_by(year) %>%
  mutate(value = na_interpolation(cases))


hosp_days<-which(time_vec%in%df1$date)


# check original data vs interpolated
#plot(df0$date, df0$cases)
#lines(df1$date,df1$value)

# ggplot(data = df1, aes(x=date,y=value))+
#  geom_bar(stat = "identity")+
#  geom_line(data=total_cases_avg,aes(x=date,y=cases))
#   

#### hosp cases series
hosp_cases<-data.frame(
  cases= df1$value,
  day  = seq(from=hosp_days[1],by=7,length.out=length(df1$date)) )


hosp_cases$date<-df1$date

id<-which(hosp_cases$date<epi_week("2020-01-01")) # fit to before covid

# Keep for plotting
hosp_cases<-hosp_cases[id,]
hosp_cases$week<-as.factor(isoweek(hosp_cases$date))
# data object
df_hosp_all<-data.frame(
  time=hosp_cases$day,
  hosp=hosp_cases$cases
)


# Bring all data together-----------------------------------------------

temp<-data.frame(time=seq(0,length(time_vec)-1))

data_long<-merge(temp,df_iid2, by.x="time", by.y = "time", all.x = TRUE)
data_long<-merge(data_long,df_sgss_all, by.x="time", by.y = "time", all.x = TRUE)
data_long<-merge(data_long,df_sgss_avg, by.x="time", by.y = "time", all.x = TRUE)
data_long<-merge(data_long,df_hosp_all, by.x="time", by.y = "time", all.x = TRUE)
data_long<-merge(data_long,df_sgss_strain, by.x="time", by.y = "time", all.x = TRUE)
data_long<-merge(data_long,df_sgss_age, by.x="time", by.y = "time", all.x = TRUE)
data_long<-merge(data_long,df_gii4_kids, by.x="time", by.y = "time", all.x = TRUE)


data_short<-data_long[which(rowSums(data_long[,2:ncol(data_long)],na.rm=TRUE)>0), ]



# Save data-----------------------------------------------------------------

qs_save(data_long, outfile1)
qs_save(data_short, outfile2)

data_for_plots<-list(
  data_iid2.c4=data_iid2.c4,
  agg_age=agg_age,
  total_cases=total_cases,
  total_cases_avg=total_cases_avg,
  hosp_cases=hosp_cases,
  total_cases_str=total_cases_str,
  dfsero=dfsero)
  
  

qs_save(data_for_plots,outfile3)


