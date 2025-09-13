# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------

root    <- here::here()
infile0  <- file.path(root,"output", "parameters.qs2")
infile  <- file.path(root, "output", "school_uk.qs2")
outfile <- file.path(root, "output", "covid_sche.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))



# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------
COVID <- "ON"

parameters <- qs_read(infile0)

tstep<-parameters$tstep
time_vec<-parameters$time_vec

school_uk<- qs2::qs_read(infile)

## Lockdown dates
covid_start<-epi_week("2020-03-23")
covid_end  <-epi_week("2021-03-16")
covid_vec<-seq(covid_start,covid_end,tstep)
id_day_covid<-which(time_vec%in%covid_vec)

covid_sche<-school_uk*0

if (COVID == "ON"){
  
  covid_sche[id_day_covid]<-1
  
  
  # Lockdown 1 = 23rd March -  3rd June 2020
  lock_start1<-epi_week("2020-03-23")
  lock_end1  <-epi_week("2020-06-03")
  lock_vec<-seq(lock_start1,lock_end1,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-1
  
  # Lockdown 1 easing =  4th June - 29th July 2020
  lock_start2<-epi_week("2020-06-04")
  lock_end2  <-epi_week("2020-07-29")
  lock_vec<-seq(lock_start2,lock_end2,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-2
  
  # Reduced restrictions = 30th July - 3rd Sep 2020
  lock_start3<-epi_week("2020-07-30")
  lock_end3  <-epi_week("2020-09-03")
  lock_vec<-seq(lock_start3,lock_end3,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-3
  
  # Schools open = 4th Sept - 26th October 2020
  lock_start4<-epi_week("2020-09-04")
  #lock_end4  <-epi_week("2020-10-26")
  lock_end4  <-epi_week("2020-11-04")
  lock_vec<-seq(lock_start4,lock_end4,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-4
  
  # Lockdown 2 = 5th November - 2nd December 2020
  
  lock_start5<-epi_week("2020-11-05") #
  lock_end5  <-epi_week("2020-12-02")
  lock_vec<-seq(lock_start5,lock_end5,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-5
  
  # Lockdown 2 easing = 3rd December - 19th December 2020
  lock_start6<-epi_week("2020-12-03") 
  lock_end6  <-epi_week("2020-12-19")
  lock_vec<-seq(lock_start6,lock_end6,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-6
  
  # Christmas = 20 December 2020 - 2nd January 2021
  lock_start7<-epi_week("2020-12-20") 
  lock_end7  <-epi_week("2021-01-04")#<- changed for continuity
  lock_vec<-seq(lock_start7,lock_end7,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-7
  
  # Lockdown 3 = 5th January - 8th March 2021
  lock_start8<-epi_week("2021-01-05") 
  lock_end8  <-epi_week("2021-03-08")
  lock_vec<-seq(lock_start8,lock_end8,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-8
  
  # Lockdown 3 with schools open = 8th March - 16th March 2021
  lock_start9<-epi_week("2021-03-09") 
  #lock_end9  <-epi_week("2021-03-16")
  lock_end9  <-epi_week("2021-06-01")
  lock_vec<-seq(lock_start9,lock_end9,tstep)
  id_day_lock<-which(time_vec%in%lock_vec)
  covid_sche[id_day_lock]<-9
}

# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(covid_sche, outfile)
