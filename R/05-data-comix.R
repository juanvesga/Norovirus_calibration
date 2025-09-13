
# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------

root    <- here::here()
infile0  <- file.path(root,"output", "parameters.qs2")
infile  <- file.path(root, "data", "raw", "contact_matrices_9_periods.csv")
infile2  <- file.path(root, "output", "polymod.qs2")
outfile <- file.path(root, "output", "comix.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))



# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------

parameters <- qs_read(infile0)

dat <- read.csv(infile, header=TRUE)
dat[,1]<-NULL
comix.period <- names(table(dat$period))
names(dat) <- c("contactee","contact","vals","period")  

p1<-t(matrix(dat[dat$period =="1. Lockdown 1",3], ncol = length(unique(dat$contactee))))
p2<-t(matrix(dat[dat$period =="2. Lockdown 1 easing",3], ncol = length(unique(dat$contactee))))
p3<-t(matrix(dat[dat$period =="3. Relaxed restrictions",3], ncol = length(unique(dat$contactee))))
p4<-t(matrix(dat[dat$period =="4. School reopening",3], ncol = length(unique(dat$contactee))))
p5<-t(matrix(dat[dat$period =="5. Lockdown 2",3], ncol = length(unique(dat$contactee))))
p6<-t(matrix(dat[dat$period =="6. Lockdown 2 easing",3], ncol = length(unique(dat$contactee))))
p7<-t(matrix(dat[dat$period =="7. Christmas",3], ncol = length(unique(dat$contactee))))
p8<-t(matrix(dat[dat$period =="8. Lockdown 3",3], ncol = length(unique(dat$contactee))))
p9<-t(matrix(dat[dat$period =="9. Lockdown 3 + schools",3], ncol = length(unique(dat$contactee))))

conmat<- qs_read(infile2)

cmx_1<-conmat$transmission*0
cmx_2<-conmat$transmission*0
cmx_3<-conmat$transmission*0
cmx_4<-conmat$transmission*0
cmx_5<-conmat$transmission*0
cmx_6<-conmat$transmission*0
cmx_7<-conmat$transmission*0
cmx_8<-conmat$transmission*0
cmx_9<-conmat$transmission*0
dat.upper.age<-c(4,11,17,29,39,49,59,69,100)
cmat.upper.age<-parameters$ages

for (ii in 1:length(parameters$ages)){
  
  ii.cmx<- min(which(dat.upper.age>=cmat.upper.age[ii]))
  
  for (jj in 1:length(parameters$ages)){
    
    jj.cmx<- min(which(dat.upper.age>=cmat.upper.age[jj]))
    
    
    
    cmx_1[ii,jj]<-p1[ii.cmx,jj.cmx]
    cmx_2[ii,jj]<-p2[ii.cmx,jj.cmx]
    cmx_3[ii,jj]<-p3[ii.cmx,jj.cmx]
    cmx_4[ii,jj]<-p4[ii.cmx,jj.cmx]
    cmx_5[ii,jj]<-p5[ii.cmx,jj.cmx]
    cmx_6[ii,jj]<-p6[ii.cmx,jj.cmx]
    cmx_7[ii,jj]<-p7[ii.cmx,jj.cmx]
    cmx_8[ii,jj]<-p8[ii.cmx,jj.cmx]
    cmx_9[ii,jj]<-p9[ii.cmx,jj.cmx]
    
    
  }
}

pop<-parameters$pop

comix<-list(
cmx_1=symm_mat(cmx_1)/rep(pop, each = ncol(conmat$transmission)),
cmx_2=symm_mat(cmx_2)/rep(pop, each = ncol(conmat$transmission)),
cmx_3=symm_mat(cmx_3)/rep(pop, each = ncol(conmat$transmission)),
cmx_4=symm_mat(cmx_4)/rep(pop, each = ncol(conmat$transmission)),
cmx_5=symm_mat(cmx_5)/rep(pop, each = ncol(conmat$transmission)),
cmx_6=symm_mat(cmx_6)/rep(pop, each = ncol(conmat$transmission)),
cmx_7=symm_mat(cmx_7)/rep(pop, each = ncol(conmat$transmission)),
cmx_8=symm_mat(cmx_8)/rep(pop, each = ncol(conmat$transmission)),
cmx_9=symm_mat(cmx_9)/rep(pop, each = ncol(conmat$transmission))
)


# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(comix, outfile)


