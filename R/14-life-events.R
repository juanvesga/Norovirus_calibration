# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")
infile_model2  <- file.path(root,"models", "model_2pars.R")
infile_chains  <- file.path(root, "output", "mcmc_chains.qs2")
infile_polymod <- file.path(root, "output", "polymod_unweighted.qs2")
infile_runs      <- file.path(root,"output", "mcmc_fits_2.qs2")

outfile1        <- file.path(root, "output", "mcmc_fits.qs2")
outfile2        <- file.path(root, "output", "mcmc_inits.qs2")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "r0_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(abind, include.only = c("abind"))


library(ggplot2)
library(dplyr)
library(matrixStats)


nsamps<- 100
step_yr<-365

# Source Odin model
parameters <- qs_read(infile0)
pars_list  <- qs_read(infile_input)
chains     <- qs_read(infile_chains) 

#combine chains
longchain<-do.call(rbind, chains)

#MLE
id<-which(longchain[,"logDensity"]==max(longchain[,"logDensity"]))[1]

exp(longchain[id,])

#get last sample
params<- exp(longchain[sample(nrow(longchain),size=nsamps,replace=FALSE),1:15])

sims     <-qs_read(infile_runs)
idx<- seq(1,dim(sims)[1],1) 
names(idx)<-paste(rownames(sims))



agg_ages<-function(index,sim){
  

  id<-c(idx[paste("inc_gi3",as.character(index[1]),sep="_")],
        idx[paste("inc_gi",as.character(index[1]),sep="_")],
        idx[paste("inc_gii4",as.character(index[1]),sep="_")],
        idx[paste("inc_gii",as.character(index[1]),sep="_")])
  
  id2<- idx[paste("pop_all",as.character(index[1]),sep="")]
  times<- (step_yr*2):dim(sim)[3]
  
  inc<-sim[id,,times]
  pop<-sim[id2,,times]
  
  
  if (length(index)>1){
    for (ii in 2:length(index)){
     
      id<-c(idx[paste("inc_gi3",as.character(index[ii]),sep="_")],
            idx[paste("inc_gi",as.character(index[ii]),sep="_")],
            idx[paste("inc_gii4",as.character(index[ii]),sep="_")],
            idx[paste("inc_gii",as.character(index[ii]),sep="_")])
      
      id2<- idx[paste("pop_all",as.character(index[1]),sep="")]
      
      inc<- inc + sim[id,,times]
      pop<- pop + sim[id2,,times] 
      
    }
  }
  
  tmp1<-colSums(inc,dims = 1 )
  
  out<-list(
    inc=t(apply(tmp1, 1, cumsum)),
    pop=pop
  )
  
  return(out)
  
}

agg_ages_strain<-function(index,sim,strain){
  
  id<-idx[paste("inc",strain,as.character(index[1]),sep="_")]
  
  id2<- idx[paste("pop_all",as.character(index[1]),sep="")]
  times<- (step_yr*2):dim(sim)[3]
  
  inc<-sim[id,,times]
  pop<-sim[id2,,times]
  
  
  if (length(index)>1){
    for (ii in 2:length(index)){
      id<-idx[paste("inc",strain,as.character(index[ii]),sep="_")]
      
      id2<- idx[paste("pop_all",as.character(index[ii]),sep="")]
      
      inc<- inc + sim[id,,times]
      pop<- pop + sim[id2,,times] 
      
    }
  }
  
  tmp1<-inc
  
  out<-list(
    inc=t(apply(tmp1, 1, cumsum)),
    pop=pop
  )
  
  return(out)
  
}

browser()

index<-1:4
inc_a0_4<-agg_ages(index,sims)$inc
inc_a0_4_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a0_4_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a0_4_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a0_4_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a0_4<-agg_ages(index,sims)$pop

le<-2 # time point 
infpc_0_4<-inc_a0_4[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gi3<-inc_a0_4_gi3[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gi<-inc_a0_4_gi[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gii4<-inc_a0_4_gii4[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gii<-inc_a0_4_gii[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])

##
index<-5:8
inc_a5_15<-agg_ages(index,sims)$inc
inc_a5_15_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a5_15_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a5_15_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a5_15_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a5_15<-agg_ages(index,sims)$pop
le<-5
infpc_5_15<-inc_a5_15[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gi3<-inc_a5_15_gi3[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gi<-inc_a5_15_gi[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gii4<-inc_a5_15_gii4[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gii<-inc_a5_15_gii[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])


index<-9:11
inc_a16_45<-agg_ages(index,sims)$inc
inc_a16_45_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a16_45_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a16_45_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a16_45_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a16_45<-agg_ages(index,sims)$pop

le<-15
infpc_16_45<-inc_a16_45[,step_yr*le]/rowMeans(pop_a16_45[,1:step_yr*le])
infpc_16_45_gi3<-inc_a16_45_gi3[,step_yr*le]/rowMeans(pop_a16_45[,1:step_yr*le])
infpc_16_45_gi<-inc_a16_45_gi[,step_yr*le]/rowMeans(pop_a16_45[,1:step_yr*le])
infpc_16_45_gii4<-inc_a16_45_gii4[,step_yr*le]/rowMeans(pop_a16_45[,1:step_yr*le])
infpc_16_45_gii<-inc_a16_45_gii[,step_yr*le]/rowMeans(pop_a16_45[,1:step_yr*le])


index<-12:13
inc_a46_65<-agg_ages(index,sims)$inc
inc_a46_65_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a46_65_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a46_65_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a46_65_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a46_65<-agg_ages(index,sims)$pop

le<-10
infpc_46_65<-inc_a46_65[,step_yr*le]/rowMeans(pop_a46_65[,1:step_yr*le])
infpc_46_65_gi3<-inc_a46_65_gi3[,step_yr*le]/rowMeans(pop_a46_65[,1:step_yr*le])
infpc_46_65_gi<-inc_a46_65_gi[,step_yr*le]/rowMeans(pop_a46_65[,1:step_yr*le])
infpc_46_65_gii4<-inc_a46_65_gii4[,step_yr*le]/rowMeans(pop_a46_65[,1:step_yr*le])
infpc_46_65_gii<-inc_a46_65_gii[,step_yr*le]/rowMeans(pop_a46_65[,1:step_yr*le])


index<-14
inc_a65p<-agg_ages(index,sims)$inc
inc_a65p_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a65p_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a65p_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a65p_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a65p<-agg_ages(index,sims)$pop

le<-7
infpc_65p<-inc_a65p[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gi3<-inc_a65p_gi3[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gi<-inc_a65p_gi[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gii4<-inc_a65p_gii4[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gii<-inc_a65p_gii[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])

len=length(infpc_5_15)
df1<-data.frame(
  y = c(infpc_0_4_gi3,
        infpc_0_4_gi3+infpc_5_15_gi3,
        infpc_0_4_gi3+infpc_5_15_gi3+infpc_16_45_gi3,
        infpc_0_4_gi3+infpc_5_15_gi3+infpc_16_45_gi3+infpc_46_65_gi3,
        infpc_0_4_gi3+infpc_5_15_gi3+infpc_16_45_gi3+infpc_46_65_gi3+infpc_65p_gi3),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 45",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("GI.3",len*5)
  
)

df2<-data.frame(
  y = c(infpc_0_4_gi,
        infpc_0_4_gi+infpc_5_15_gi,
        infpc_0_4_gi+infpc_5_15_gi+infpc_16_45_gi,
        infpc_0_4_gi+infpc_5_15_gi+infpc_16_45_gi+infpc_46_65_gi,
        infpc_0_4_gi+infpc_5_15_gi+infpc_16_45_gi+infpc_46_65_gi+infpc_65p_gi),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 45",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("Other GI",len*5)
  
)

df3<-data.frame(
  y = c(infpc_0_4_gii4,
        infpc_0_4_gii4+infpc_5_15_gii4,
        infpc_0_4_gii4+infpc_5_15_gii4+infpc_16_45_gii4,
        infpc_0_4_gii4+infpc_5_15_gii4+infpc_16_45_gii4+infpc_46_65_gii4,
        infpc_0_4_gii4+infpc_5_15_gii4+infpc_16_45_gii4+infpc_46_65_gii4+infpc_65p_gii4),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 45",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("GII.4",len*5)
  
)

df4<-data.frame(
  y = c(infpc_0_4_gii,
        infpc_0_4_gii+infpc_5_15_gii,
        infpc_0_4_gii+infpc_5_15_gii+infpc_16_45_gii,
        infpc_0_4_gii+infpc_5_15_gii+infpc_16_45_gii+infpc_46_65_gii,
        infpc_0_4_gii+infpc_5_15_gii+infpc_16_45_gii+infpc_46_65_gii+infpc_65p_gii),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 45",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("Other GII",len*5)
  
)

df5<-data.frame(
  y = c(infpc_0_4,
        infpc_0_4+infpc_5_15,
        infpc_0_4+infpc_5_15+infpc_16_45,
        infpc_0_4+infpc_5_15+infpc_16_45+infpc_46_65,
        infpc_0_4+infpc_5_15+infpc_16_45+infpc_46_65+infpc_65p),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 45",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("All",len*5)
  
)


df<-rbind(df1,df2,df3,df4,df5)
df$Age<-factor(df$Age,
               levels=c("at 4","at 15","at 45","at 65","at 80"))
df$Strain<-factor(df$Strain,
                  levels=c("GI.3","Other GI","GII.4","Other GII","All"))

library(dplyr)
library(viridis)
plot_strains<-df %>%
  ggplot( aes(x=Age, y=y, fill=Strain)) +
  #geom_jitter(aes(color=Strain), size=0.4, alpha=0.5) +
  geom_boxplot(outlier.shape = NA)  +
  geom_boxplot(aes(color = Strain),
               fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
               show.legend = F)+
  theme_light() +
  scale_fill_manual(values = c('skyblue3','firebrick2','yellow3','grey40','slateblue2'))+
  scale_color_manual(values = c('skyblue3','firebrick2','yellow3','grey40','slateblue2'))+
  #ggtitle("Accrued Norovirus infection") +
  ylim(0,4)+
  xlab("Age") + ylab("AGE episodes")  +
  theme(
    #legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10), legend.key = element_blank(),
    # plot.tag = element_text(size = 12, face = "bold"),
    #  plot.tag.position = c(0.1, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )


windows()
plot_strains


len=length(infpc_5_15)
df1<- data.frame(y=c(infpc_0_4,
                     infpc_0_4+infpc_5_15, 
                     infpc_0_4+infpc_5_15+infpc_16_45,
                     infpc_0_4+infpc_5_15+infpc_16_45+infpc_46_65,  
                     infpc_0_4+infpc_5_15+infpc_16_45+infpc_46_65+ infpc_65p),
                 age=c(rep("at 4",len),
                       rep("at 15",len),
                       rep("at 45",len),
                       rep("at 65",len),
                       rep("at 80",len)))

df1$age<-factor(df1$age,
                levels=c("at 4","at 15","at 45","at 65","at 80")) 


plot_all<-df1 %>%
  ggplot( aes(x=age, y=y, fill=age)) +
  geom_jitter(color="grey", size=0.4, alpha=0.5) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_minimal() +
  ggtitle("Accrued Norovirus infection") +
  xlab("Age") + ylab("Mean number of symptomatic events")                 

windows()
gridExtra::grid.arrange(plot_strains)
windows()
gridExtra::grid.arrange(plot_all)

rm(sims,out)

