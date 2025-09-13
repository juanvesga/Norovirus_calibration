
# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------

root    <- here::here()
infile0  <- file.path(root,"output", "parameters.qs2")
infile  <- file.path(root, "data", "raw", "uk_holidays.csv")
outfile <- file.path(root, "output", "school_uk.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(lubridate, include.only = c("ceiling_date","leap_year"))

# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------

parameters <- qs_read(infile0)

dat <- read.csv(infile, header=TRUE)

id<-which(dat$Date=="28-Feb")

school_year<-dat$School

val<-school_year[id]

school_year_leap <-c(school_year[1:id], val, school_year[- (1:id)])

projection_enddate  <-epi_week("2035-12-28")#epi_week("2023-12-31")
sche_vec<-seq(parameters$sim_startdate,projection_enddate,"year")
isleap<-leap_year(sche_vec)

sche_list<-list()

for (gg in 1:length(sche_vec)){
  
  if (isleap[gg]){
    sche_list[[gg]]<-school_year_leap
  }else{
    sche_list[[gg]]<-school_year
  }
}

school_uk<-unlist(sche_list)



# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(school_uk, outfile)






