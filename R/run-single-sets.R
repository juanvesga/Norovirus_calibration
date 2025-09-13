# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_dataS   <- file.path(root,"output", "data_short.qs2")
infile_dataL   <- file.path(root,"output", "data_long.qs2")
infile_dataPlots <- file.path(root, "output", "data_for_plots.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")
infile_prior   <- file.path(root,"output", "priors.qs2")
infile_samples <-  file.path(root, "output", "chains.qs")
infile_index <-  file.path(root, "output", "index.qs")
infile_model2  <- file.path(root,"models", "model_2pars.R")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "plot_functions.R"))


# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(odin2)
library(dust2)
library(monty)
library(posterior)
library(bayesplot)
library(matrixStats)
library(ggplot2)
#library(introdataviz)
library(see)
library(tidyr)

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------
index    <-qs_read(infile_index)
data_all <-qs_read(infile_dataS)
data_plots<-qs_read(infile_dataPlots)

tt<-data_all$time
total_cases_str<-data_plots$total_cases_str
total_cases<-data_plots$total_cases
agg_age<-data_plots$agg_age
data_iid2.c4<-data_plots$data_iid2.c4
dfsero<-data_plots$dfsero

# Source Odin model
source(infile_model2)
parameters <- qs_read(infile0)
data       <- qs_read(infile_dataS)
pars_list  <- qs_read(infile_input)
prior      <- qs_read(infile_prior)
#Parameters object
pars<-pars_list$pars_list
packer<-pars_list$packer
start_list<-pars_list$start_pars



#dust object
filter <- dust_unfilter_create(noro_model(), 
                                   time_start = 0, 
                                   data = data,
                               n_particles = 2)


## Choose parameters
# start_pars<-list(
#   beta_1      =    0.13410824,
#   beta_2      =    0.12708480,
#   beta_3      =    0.20300807,
#   beta_4      =    0.27545582,
#   aduRR       =    0.58903565,
#   maternalAB  =  367.16148967,
#   imm_yr      =     16.15967444 ,
#   imm_fac     =     1.70697752,
#   repfac_0    =     64.26091774,
#   repfac_5    =     450.23810289,
#   repfac_15   =     345.14095794,
#   repfac_65p  =    59.89184305,
#   crossp_GI   =        0.09838119,
#   crossp_GII  =      0.01773758 )

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
crossp_GI=0.04103662   ,
crossp_GII=0.01982703
)

# start_pars<-list(
#   beta_1    = 0.13953720,   
#   beta_2    = 0.12680114,          
#   beta_3    = 0.15410557,          
#   beta_4    = 0.21306245,           
#   aduRR     = 0.369372532  , 
#   maternalAB= 160.61197156,  
#   imm_yr    = 5.18686356,   
#   imm_fac   = 1.68360017 ,
#   repfac_0  = 105.07189709,  
#   repfac_5  = 954.17803464 ,    
#   repfac_15 = 581.94243972 , 
#   repfac_65p= 95.41216406  , 
#   crossp_GI = 0.04103662   ,
#   crossp_GII= 0.01982703   
# )



pars2=c(start_pars,pars_list$fixed_pars)

dust_likelihood_run(filter, pars2, save_trajectories = TRUE)

state <- dust_likelihood_last_trajectories(filter)
state <- aperm(state, c(1,3,2))



# SGSS reported time series -----------------------------------------------



ii<-which(!is.na(data_all$reported_gi3))

# GI3
dat <- data.frame(
  x = total_cases_str$x,
  y=data_all$reported_gi3[ii]
)

fitgi3<-plot_reported_strain(state[index$reported_wk_gi3,ii,],dat,1)

dat <- data.frame(
  x = total_cases_str$x,
  y=data_all$reported_gi[ii]
)
fitgi<-plot_reported_strain(state[index$reported_wk_gi,ii,],dat,2)

dat <- data.frame(
  x = total_cases_str$x,
  y=data_all$reported_gii4[ii]
)

fitgii4<-plot_reported_strain(state[index$reported_wk_gii4,ii,],dat,3)

dat <- data.frame(
  x = total_cases_str$x,
  y=data_all$reported_gii4[ii]
)
fitgii<-plot_reported_strain(state[index$reported_wk_gii,ii,],dat,4)




reported<-t(colSums(state[index$reported_wk, , ]))

ii<-which(!is.na(data_all$reported))


df<-data.frame(t(reported[,ii]))

df$x<-total_cases$date 
dat_sim <- reshape2::melt(df, id = "x")

df_d <- data.frame(
  x = total_cases$date,
  y=data_all$reported[ii]
)

df_s <- as.data.frame(
  rowQuantiles(t(reported[,ii]),
               probs = c(0.025, 0.5, 0.975)
  )
)
df_s$x<-total_cases$date

col_0<-"#cc0044"
col_1<-"#E54C20"
col_2<-"#8c8cd9"
col_3<-"#00A04B"
col_4<-"grey40"

start_date <- total_cases$date[1]# ymd("2014-01-01")
end_date <- total_cases$date[length(total_cases$date)]# ymd("2020-31-12")

fits_sgss <- ggplot() +
  geom_line(data = dat_sim, aes(x = x, y=value, group=variable),
            col = col_3 , alpha = 0.1, lwd = 0.001) +
  geom_line(data = df_s, aes(x = x, y=`50%`), col = "black",lwd = 1) +
  geom_bar(data = df_d, aes(x = x, y = y), stat="identity", #width = 0.9, 
           fill=col_2, alpha=0.5) +
  labs(tag = "A", x = "", y = "Weekly \n cases reported") +
  theme_classic() +
  ylim(c(0,max(df_d$y)*1.5))+
  scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y", 
               limits = c( start_date, end_date), expand=c(0,0))+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    axis.text.x = element_text(size=8,angle = 60, hjust = 1),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.1, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )

#}




# IID incidence by strain and age -----------------------------------------

nice_cols = c( "#007A87","#FF5A5F","#FFB400", "#7B0051", 
               "#8CE071",  "#00D1C1", "#FFAA91", "#B4A76C", 
               "#9CA299", "#565A5C", "#00A04B", "#E54C20")



ii<-which(!is.na(data_all$cases_a1))


modelled_irate <-rbind(
  ( (    state[index$inc_year_gi3[1],ii,]+
           state[index$inc_year_gi[1],ii,]+
           state[index$inc_year_gii4[1],ii,]+
           state[index$inc_year_gii[1],ii,])/(state[index$pop_by4age[1],ii,])),
  
  ( (    state[index$inc_year_gi3[2],ii,]+
           state[index$inc_year_gi[2],ii,]+
           state[index$inc_year_gii4[2],ii,]+
           state[index$inc_year_gii[2],ii,])/(state[index$pop_by4age[2],ii,])),
  
  ( (    state[index$inc_year_gi3[3],ii,]+
           state[index$inc_year_gi[3],ii,]+
           state[index$inc_year_gii4[3],ii,]+
           state[index$inc_year_gii[3],ii,])/(state[index$pop_by4age[3],ii,])),
  
  ( (    state[index$inc_year_gi3[4],ii,]+
           state[index$inc_year_gi[4],ii,]+
           state[index$inc_year_gii4[4],ii,]+
           state[index$inc_year_gii[4],ii,])/(state[index$pop_by4age[4],ii,])),
  
  ( (    state[index$inc_year_gi3[5],ii,]+
           state[index$inc_year_gi[5],ii,]+
           state[index$inc_year_gii4[5],ii,]+
           state[index$inc_year_gii[5],ii,])/(state[index$pop_by4age[5],ii,]))
)*1000



iid2_all<-iid2_plot_func(modelled_irate, data_iid2.c4,   "all"," ", "F",nice_cols[8])




# GII4 prevalence in children  --------------------------------------------


x_d <- c(1, 2, 3, 4, 5, 6)-0.28 

ii<-which(!is.na(data_all$sero1))

observed_size<-c(
  103,
  107,
  121,
  124,
  122,
  109
)
sero_obs<-(c(data_all$sero1[ii], data_all$sero2[ii], data_all$sero3[ii], data_all$sero4[ii],
             data_all$sero5[ii], data_all$sero6[ii])/observed_size)*1e2

df_d <- data.frame(
  x<-x_d,
  sero=dfsero$mean*100,
  low=dfsero$low*100,
  up=dfsero$up*100)

fac<-c(2,3,4,5,6,7)


#id<-which(tt%in%which(!is.na(data$sero1)))

sero_model<-state[index$seroprev[2:7], ,ii]*100



df_qtls <- as.data.frame(rowQuantiles((sero_model),
                                      probs = c(0.025, 0.5, 0.975)))

df1 <- data.frame(t(sero_model)) 
colnames(df1) <- paste(c("[0 1)","[1 2)","[2 3)","[3 4)","[4 5)","[5 6)"))
df_m <- reshape2::melt(df1)
df_m$variable <- as.factor(df_m$variable)
#df_d$x <- factor(c("[0 1)","[1 2)","[2 3)","[3 4)","[4 5)","[5 6)"))
df_qtls$x <- factor(c("[0 1)","[1 2)","[2 3)","[3 4)","[4 5)","[5 6)"))


viol_col <-  nice_cols[3]
err_col <- "black"
data_col <- "black"



fits_sero <-  ggplot() +
  geom_point(
    data = df_m,
    aes(x = variable, y = value),
    color=viol_col,
    position = position_jitter(w = .15), 
    size = 1,
    alpha = 0.15,
    show.legend = F) +
  geom_point(data = df_d, mapping = aes(x = x, y = sero, color = "Data (95% CI)"), size = 2, shape = 15) +
  geom_errorbar(
    mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
    width = .15, position = position_dodge(.5)
  ) +
  geom_boxplot(
    data = df_m,
    aes(x = variable, y = value, fill = "Posterior Density"),
    width=0.15,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  labs(tag = "H", 
       x = "Age", 
       y = "GII.4 Seroprevalence (%)") +
  theme_classic() +
  
  ylim(0, max(df_d$up)*1.2) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    # axis.text.x = element_text(angle = 60, hjust = 1),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.2, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )




# SGSS age reporting fractions --------------------------------------------


ii<-which(!is.na(data_all$a0_event))


a0_model<-((state[index$infections_day_gi3[1],ii,]+
              state[index$infections_day_gi[1],ii,]+
              state[index$infections_day_gii4[1],ii,]+
              state[index$infections_day_gii[1],ii,]) /
             (    
               colSums(state[index$infections_day_gi3,ii,])+
                 colSums(state[index$infections_day_gi,ii,])+
                 colSums(state[index$infections_day_gii4,ii,])+
                 colSums(state[index$infections_day_gii,ii,])))


a5_model<-((state[index$infections_day_gi3[2],ii,]+
              state[index$infections_day_gi[2],ii,]+
              state[index$infections_day_gii4[2],ii,]+
              state[index$infections_day_gii[2],ii,]) /
             (    
               colSums(state[index$infections_day_gi3,ii,])+
                 colSums(state[index$infections_day_gi,ii,])+
                 colSums(state[index$infections_day_gii4,ii,])+
                 colSums(state[index$infections_day_gii,ii,])))

a15_model<-((state[index$infections_day_gi3[3],ii,]+
               state[index$infections_day_gi[3],ii,]+
               state[index$infections_day_gii4[3],ii,]+
               state[index$infections_day_gii[3],ii,]) /
              (    
                colSums(state[index$infections_day_gi3,ii,])+
                  colSums(state[index$infections_day_gi,ii,])+
                  colSums(state[index$infections_day_gii4,ii,])+
                  colSums(state[index$infections_day_gii,ii,])))

a65_model<-((state[index$infections_day_gi3[4],ii,]+
               state[index$infections_day_gi[4],ii,]+
               state[index$infections_day_gii4[4],ii,]+
               state[index$infections_day_gii[4],ii,]) /
              (    
                colSums(state[index$infections_day_gi3,ii,])+
                  colSums(state[index$infections_day_gi,ii,])+
                  colSums(state[index$infections_day_gii4,ii,])+
                  colSums(state[index$infections_day_gii,ii,])))





age_model<-cbind(a0_model,a5_model, a15_model,a65_model )*100


x_d <- c(1,2,3,4)-0.28 # bin x axis positions


df_d <- data.frame(
  x = x_d,
  age = agg_age[,2]*100,
  low = agg_age[,1]*100,
  up = agg_age[,3]*100
)

df_qtls <- as.data.frame(rowQuantiles(t(age_model),
                                      probs = c(0.025, 0.5, 0.975)))

df1 <- data.frame((age_model)) 
colnames(df1) <- paste(c("[0 4)","[5 14)","[15 64)", "65+"))
df_m <- reshape2::melt(df1)
df_m$variable <- as.factor(df_m$variable)
#df_d$x <- factor(c("[0 4)","[5 14)","[15 64)", "65+"))
df_qtls$x <- factor(c("[0 4)","[5 14)","[15 64)", "65+"))


viol_col <-  nice_cols[9]
err_col <- "black"
data_col <- "black"

fits_age <- ggplot() +
  geom_point(
    data = df_m,
    aes(x = variable, y = value),
    color=viol_col,
    position = position_jitter(w = .15), 
    size = 1,
    alpha = 0.15,
    show.legend = F) +
  geom_point(data = df_d, 
             mapping = aes(x = x, y = age, color = "Data (95% CI)"),
             size = 2, shape = 15) +
  geom_errorbar(
    mapping = aes(x = x, ymin = low, ymax = up), 
    data = df_d,
    width = .15, position = position_dodge(.5)
  ) +
  # geom_split_violin(
  #   data = df_m,
  #   aes(x = variable, y = value, fill = "Posterior Density"),
  #   alpha = .4, 
  #   trim = FALSE
  # ) +
  geom_boxplot(
    data = df_m,
    aes(x = variable, y = value, fill = "Posterior Density"),
    width=0.15,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  labs(tag = "G", x = "Age", y = "Reported (%)") +
  theme_classic() +
  
  ylim(0, max(df_d$up)*1.2) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    # axis.text.x = element_text(angle = 60, hjust = 1),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.2, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )





# Get plots together  -----------------------------------------------------

windows()


gridExtra::grid.arrange(
  fits_sgss, 
  fitgi3,
  fitgi,
  fitgii4,
  fitgii,
  iid2_all,
  fits_age, 
  fits_sero,
  layout_matrix = rbind(c(1 ,1),
                        c(2, 3),
                        c(4, 5),
                        c(6, 7),
                        c(8, 8)))




