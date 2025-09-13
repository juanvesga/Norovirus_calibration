# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile0  <- file.path(root,"output", "parameters.qs2")
infile  <- file.path(root, "data", "raw", "vanHoek.csv")
outfile <- file.path(root, "output", "polymod.qs2")
outfile2 <- file.path(root, "output", "polymod_unweighted.qs2")

# Packages
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "modify_attach.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(socialmixr, include.only = c("contact_matrix"))



# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------
ages   <-  c(1,2,3,4,5,6,7,15,25,35,45,55,65,75) # Upper end of age bands

age_sq <- ages-1

adults <-  seq(tail(which(ages<=5),1),length(ages),1) 

infa_id<- which(ages<5)

adult_id<-which(ages>=15)

da <- diff(c(0,ages))

aging_vec = c(1/head(da,-1),0)

aging_vec[]

# find diagonal matrix for aging transitions
aging_mat <- diag(-1/da)
aging_mat[row(aging_mat)-col(aging_mat)==1] <- 1/head(da,-1)
aging_mat[length(ages),length(ages)]<-0 # Last age group aging accounted for with mortality
age.categories <- as.factor(ages)

#Load age contact matrix for Infants in UK according to vanHoek supplement

data(polymod, package = "socialmixr") # POLYMOD for all other contacts

contact = contact_matrix(
  polymod, countries = "United Kingdom", 
  age.limits = c(0,as.numeric(as.character(age.categories[1:(length(ages)-1)]))),
  symmetric = TRUE)

contact_holi = contact_matrix(
  polymod, countries = "United Kingdom", 
  age.limits = c(0,as.numeric(as.character(age.categories[1:(length(ages)-1)]))),
  filter = list(cnt_school=0),
  symmetric = TRUE)

pop = contact$demography$population # Population numbers



# Replace with vanHoek infants matrix

dat <- read.csv(infile)

vhoek<-contact$matrix*0

vhoek.upper.age<-c(1, 4,9,14,19,24,29,34,39,44,49,54,59,64,69,100)

cmat.upper.age<-ages

for (ii in 1:length(ages)){
  
  ii.cmx<- min(which( vhoek.upper.age>=cmat.upper.age[ii]))
  
  for (jj in 1:length(ages)){
    
    jj.cmx<- min(which( vhoek.upper.age>=cmat.upper.age[jj]))
    
    vhoek[ii,jj]<-dat[ii.cmx,jj.cmx]
  }
}


# Make matrix symmetric
x1<-contact$matrix
x1[1,]<-vhoek[1,]
x1[,1]<-vhoek[,1]
x<-((x1+t(x1))/2)
contact$matrix<-x

x2<-contact_holi$matrix
x2[1,]<-vhoek[1,]
x2[,1]<-vhoek[,1]
x<-((x2+t(x2))/2)
contact_holi$matrix<-x


# Matrix to input into transmission formula, note is corrected for pop size
transmission <-symm_mat(contact$matrix)/rep(pop, each = ncol(contact$matrix))

transmission_holi <- symm_mat(contact_holi$matrix)/rep(pop, each = ncol(contact$matrix))

polymod<-list()

polymod$transmission<-transmission

polymod$transmission_holi<-transmission_holi


# Uncorrected matrices

M<-list(
  m     =symm_mat(contact$matrix),
  m_holi=symm_mat(contact_holi$matrix)
)


# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
parameters <- qs_read(infile0)
parameters$pop<-pop

qs_save(polymod, outfile)
qs_save(parameters, infile0)
qs_save(M, outfile2)



