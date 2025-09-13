# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root             <- here::here()
infile0          <- file.path(root,"output", "parameters.qs2")
infile_dataS     <- file.path(root,"output", "data_short.qs2")
infile_dataL     <- file.path(root,"output", "data_long.qs2")
infile_dataPlots <- file.path(root, "output", "data_for_plots.qs2")
infile_input     <- file.path(root,"output", "params_list.qs2")
#infile_model2  <- file.path(root,"models", "model_2pars.R")
infile_model2  <- file.path(root,"models", "model_simple.R")
#infile_model2  <- file.path(root,"models", "model_simple_no_reinf.R")
#infile_model2  <- file.path(root,"models", "model_drop.R")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "plot_functions.R"))
source(file.path(root, "R", "collect_function.R"))
source(file.path(root, "R", "update_interventions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(tidyr, include.only = c("gather"))

library(odin2)
library(dust2)
library(matrixStats)
library(ggplot2)
library(see)



# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------
source(infile_model2)
data_all <-qs_read(infile_dataS)
data_plots<-qs_read(infile_dataPlots)
parameters <- qs_read(infile0)
pars_list  <- qs_read(infile_input)

# params<-c(
#   beta_1      =   0.088,# 0.088,#0.06824650,   
#   beta_2      =   0.088,# 0.0875,#0.06395356,  
#   beta_3      =   0.088*2,#,# 0.2,#0.16260249,   
#   beta_4      =   0.088,# 0.15,#0.05,#0.20521798,   
#   aduRR       =    0.1,#0.53036016, 
#   maternalAB  =    180,#388.23374771,  
#   imm_yr      =    35,#73.48567995,   
#   imm_fac     =    1.24225551, 
#   repfac_0    =    518.28692101,
#   repfac_5    =    865.87082315, 
#   repfac_15   =    629.22414099,    
#   repfac_65p  =    60.06307963,   
#   crossp_GI   =    0,#0.16961727,   
#   crossp_GII  =    0,#0.1,#0.02969143,  
#   reported_var=    0.27301538)


params<-c(
  beta_1      =    11.0117229e-02, 
  beta_2      =    11.03170465e-02, 
  beta_3      =    2.316278e-01, 
  beta_4      =    12.6375149e-02, 
  aduRR       =    0.1, 
  maternalAB  =    100, 
  imm_yr      =    60,
 # imm_fac     =   1.04,
  repfac_0    =    2.303229e+02,
  repfac_5    =    1.454593e+03, 
  repfac_15   =    6.573526e+02, 
  repfac_65p  =    3.952568e+01, 
 # reported_var=    4.889772e-01, 
  season_lag  =    1.110333e+01, 
  season_amp  =    1.052237e-01)


ref<-c(10,5)
1-(ref[1]/ sum(ref))
1-(ref[2]/ sum(ref))

1-(1-(ref[1]/ sum(ref)))


pars<-update_parameters(params,pars_list$pars_list)  

#dust object
#sys  <- dust_system_create(noro_model(), pars, n_particles = 1)

# Determinsitic object
sys  <- dust2::dust_system_create(noro_model, pars, deterministic = TRUE, n_particles = 1)

# Get state index
index<- dust2::dust_unpack_index(sys)

#Get default state initial
state<-dust2::dust_system_state(sys)

named_state<-dust2::dust_unpack_state(sys,state)

# Set inital model state based on calibrated model

# Set default initial state of the model 
dust2::dust_system_set_state_initial(sys)

# Set parameters
dust2::dust_system_update_pars(sys, pars=pars)


# run simulation
endsim<-24109
tt<-seq(0,endsim)
y0<- dust2::dust_system_simulate(sys, tt)



# unpack results
y0 <- dust2::dust_unpack_state(sys, y0)



state<-t(collect_function(y0))





# 
# 
#  matplot(t(y0$S),type="l")
# 
# 
#  plot(colSums(y0$pop_by4age),type="l")
# 
#  plot(colSums(y0$inc_day_gi3),type="l")
#  
#  plot(colSums(y0$inc_day_gi),type="l")
#  
#  plot(colSums(y0$inc_day_gii4),type="l")
#  
#  plot(colSums(y0$inc_day_gii),type="l")
#  
# 
#  plot(colSums(y0$inc_day_gi3)[1000:endsim],type="l")
# 
#  plot(colSums(y0$inc_day_gi)[1000:endsim],type="l")
# 
#  plot(colSums(y0$inc_day_gii4)[1000:endsim],type="l")
# 
# plot(colSums(y0$inc_day_gii)[1000:endsim],type="l")
# 






# Plotting functions ------------------------------------------------------


## errorbars and boxplots for IID2 incodence by age
iid2_plot_func<-function(irates, data,strain_string,title_string,tag_text,color_viol){
  
  df_qtls <- as.data.frame(rowQuantiles((irates),
                                        probs = c(0.025, 0.5, 0.975)))
  x_d <- c(1, 2, 3, 4, 5)-0.28 # bin x axis positions
  
  if (strain_string=="all"){
    df_d <- data.frame(x = x_d)
    df_d$inc = data$per1000personyears
    df_d$low = data$CI_lower
    df_d$up = data$CI_upper
  }else{
    df_d <- data.frame(x = x_d)
    df_d$inc = data[[paste0("per1000personyears","_",strain_string,sep="")]]
    df_d$low = data[[paste0("CI_lower","_",strain_string,sep="")]]
    df_d$up = data[[paste0("CI_upper","_",strain_string,sep="")]] 
  }
  
  df1 <- data.frame(t(irates)) 
  colnames(df1) <- paste(c("[0 1)","[1 4)","[5 14)","[15 64)","65+"))
  df_m <- reshape2::melt(df1)
  df_m$variable <- as.factor(df_m$variable)
  df_qtls$x <- factor(c("[0 1)","[1 4)","[5 14)","[15 64)","65+"))
  
  
  viol_col <- color_viol 
  err_col <- "black"
  data_col <- "black"
  
  iid2_plot <- ggplot() +
    geom_violin(
      data = df_m,
      aes(x = variable, y = value, fill = "Posterior Density"),
      draw_quantiles = c(0.5),
      width = 1,
      linetype = 1,
      trim = FALSE,
      color = "white",
      alpha = 0.7
    ) +
    geom_point(data = df_d, mapping = aes(x = x, y = inc, color = "Data (95% CI)"), size = 2, shape = 15) +
    geom_errorbar(
      mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
      width = .15, position = position_dodge(.5)
    ) +
    geom_boxplot(
      data = df_m,
      aes(x = variable, y = value, fill = "Posterior Density"),
      width=0.15,
      fatten=NULL,
      outlier.shape = NA,
      alpha = 0.5
    ) +
    
    labs(tag=tag_text, 
         x = "", 
         y =paste(title_string , "Incidence per 1k\n person-year",sep=" ")) +
    theme_classic() +
    
    ylim(0, max(df_d$up)*1.2) +
    scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
    scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
    theme(
      legend.position = "none",
      panel.background = element_blank(),
      axis.text = element_text(colour = "black", size = 10),
      axis.title = element_text(size = 10),
      plot.tag = element_text(size = 10),
      plot.tag.position = c(0.3, 0.91),
      panel.border = element_rect(colour = "black", fill=NA)
    )
  
  return(iid2_plot)
  
}


# SGSS reported time series -----------------------------------------------


#id<-sample(n_out*length(chain_selection),n_out)

# if (sgss_mode== "strain"){

#ii<-which(!is.na(data_all$reported_gi3))
ii<-data_all$time[which(!is.na(data_all$reported_gi3))]-1 # model strats from 0 and to match by week output

dat<-colSums(data.frame(
  GI3=data_plots$total_cases_str$gi3,
  GI=data_plots$total_cases_str$gi,
  GII4=data_plots$total_cases_str$gii4,
  GII=data_plots$total_cases_str$gii
))

dat<-data.frame(
  strain=factor(c("GI3","GI","GII4","GII"),levels=c("GI3","GI","GII4","GII")),
  y     =unlist(dat)
)


mod<-rbind(
  data.frame(y=sum(state["reported_wk_gi3",ii]),strain="GI3"),
  data.frame(y=sum(state["reported_wk_gi",ii]),strain="GI"),
  data.frame(y=sum(state["reported_wk_gii4",ii]),strain="GII4"),
  data.frame(y=sum(state["reported_wk_gii",ii]),strain="GII")
)


# plot(state["reported_wk_gi3",ii],type="l",col="blue",ylim = c(0,40))
# lines(state["reported_wk_gi",ii],col="red")
# lines(state["reported_wk_gii4",ii], col="gold3")
# lines(state["reported_wk_gii",ii], col="grey23")
# 
# 
# 
# plot(state["inc_year_gi_4",ii],type="l",col="blue",ylim=c(0,1e5))
# lines(state["inc_year_gi3_4",ii],col="red")
# lines(state["inc_year_gii4_4",ii], col="gold3")
# lines(state["inc_year_gii_4",ii], col="grey23")

sum(state["inc_year_gi3_4",ii])

sum(state["inc_year_gii4_4",ii])

mod$strain<-factor(mod$strain,levels=c("GI3","GI","GII4","GII"))


fits_bystrain <-  ggplot() +
  geom_boxplot(
    data = mod,
    aes(x = strain, y = y, fill = strain),
    width=0.15,
    #outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_point(data = dat, 
             aes(x = strain, y = y),
             position = position_dodge(0.5),
             color = "black",
             size = 3, 
             shape = 18,
             alpha =0.2) +
  
  labs(tag = "B", 
       x = "Strain", 
       y = "Cases reported") +
  theme_classic() +
  
  #ylim(0, 1) +
  scale_fill_manual(name = "", values = cols<-c("GI3"="dodgerblue",
                                                "GI" ="#cc0044",
                                                "GII4"="#FFB400",
                                                "GII"="#00A04B")) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10),
    plot.tag = element_text(size = 10),
    plot.tag.position = c(0.15, 0.91),
    panel.border = element_rect(colour = "black", fill=NA)
  )



# SGSS long series


reported<-(colSums(state[c("reported_wk_1","reported_wk_2","reported_wk_3","reported_wk_4",
                           "reported_wk_5","reported_wk_6","reported_wk_7","reported_wk_8",
                           "reported_wk_9","reported_wk_10","reported_wk_11","reported_wk_12",
                           "reported_wk_13","reported_wk_14"), ]))

ii<-which(!is.na(data_all$reported))
model_t<-data_all$time[which(!is.na(data_all$reported))]-1 # model strats from 0 and to match by week output


#df<-data.frame(t(reported[id,ii]))
df<-data.frame((reported[model_t]))

df$x<-data_plots$total_cases$date 
dat_sim <- reshape2::melt(df, id = "x")

df_d <- data.frame(
  x = data_plots$total_cases$date,
  y=data_all$reported[ii]
)


col_0<-"#cc0044"
col_1<-"#E54C20"
col_2<-"#8c8cd9"
col_3<-"#00A04B"
col_4<-"grey40"

start_date <- data_plots$total_cases$date[1]# ymd("2014-01-01")
end_date <- data_plots$total_cases$date[length(data_plots$total_cases$date)]# ymd("2020-31-12")

fits_sgss <- ggplot() +
  geom_line(data = dat_sim, aes(x = x, y=value, group=variable),
            col = col_3 , alpha = 1, lwd = 0.1) +
  geom_bar(data = df_d, aes(x = x, y = y), stat="identity", #width = 0.9, 
           fill=col_2, alpha=0.5) +
  labs(tag = "A", x = "", y = "Weekly cases\n reported") +
  theme_classic() +
  ylim(c(0,max(df_d$y)*1.5))+
  scale_x_date(date_breaks = "4 months", date_labels = "%b-%Y", 
               limits = c( start_date, end_date), expand=c(0,0))+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10),
    #legend.text = element_text(size = 9), legend.key = element_blank(),
    axis.text.x = element_text(size=8,angle = 50, hjust = 1),
    plot.tag = element_text(size = 10),
    plot.tag.position = c(0.15, 0.91),
    panel.border = element_rect(colour = "black", fill=NA)
  )

#}


# SGSS average season



ii<-which(!is.na(data_all$reported_avg))
model_t<-data_all$time[which(!is.na(data_all$reported_avg))]-1 # model strats from 0 and to match by week output


#df<-data.frame(t(reported[id,ii]))
df<-data.frame((reported[model_t]))

df$x<-data_plots$total_cases_avg$date 
dat_sim <- reshape2::melt(df, id = "x")

df_d <- data.frame(
  x = data_plots$total_cases_avg$date,
  y=data_all$reported[ii]
)


col_0<-"#cc0044"
col_1<-"#E54C20"
col_2<-"#8c8cd9"
col_3<-"#00A04B"
col_4<-"grey40"

start_date <- data_plots$total_cases_avg$date[1]# ymd("2014-01-01")
end_date <- data_plots$total_cases_avg$date[length(data_plots$total_cases_avg$date)]# ymd("2020-31-12")

fits_sgss_avg <- ggplot() +
  geom_line(data = dat_sim, aes(x = x, y=value, group=variable),
            col = col_3 , alpha = 1, lwd = 0.1) +
  geom_bar(data = df_d, aes(x = x, y = y), stat="identity", #width = 0.9, 
           fill=col_2, alpha=0.5) +
  labs(tag = "A", x = "", y = "Weekly cases\n reported") +
  theme_classic() +
  ylim(c(0,max(df_d$y)*1.5))+
  scale_x_date(date_breaks = "1 months", date_labels = "%b-%Y", 
               limits = c( start_date, end_date), expand=c(0,0))+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10),
    #legend.text = element_text(size = 9), legend.key = element_blank(),
    axis.text.x = element_text(size=8,angle = 50, hjust = 1),
    plot.tag = element_text(size = 10),
    plot.tag.position = c(0.15, 0.91),
    panel.border = element_rect(colour = "black", fill=NA)
  )




# IID incidence by strain and age -----------------------------------------

nice_cols = c( "#007A87","#FF5A5F","#FFB400", "#7B0051", 
               "#8CE071",  "#00D1C1", "#FFAA91", "#B4A76C", 
               "#9CA299", "#565A5C", "#00A04B", "#E54C20")




#ii<-which(!is.na(data_all$cases_a1))
ii<-data_all$time[which(!is.na(data_all$cases_a1))]-1 # model strats from 0 and to match by week output


modelled_irate <-rbind(
  ( (  state['inc_year_gi3_1',ii]+
         state['inc_year_gi_1',ii]+
         state['inc_year_gii4_1',ii]+
         state['inc_year_gii_1',ii])/(state['pop_by4age1',ii])),
  
  ( (  state['inc_year_gi3_2',ii]+
         state['inc_year_gi_2',ii]+
         state['inc_year_gii4_2',ii]+
         state['inc_year_gii_2',ii])/(state['pop_by4age2',ii])) ,
  
  ( (  state['inc_year_gi3_3',ii]+
         state['inc_year_gi_3',ii]+
         state['inc_year_gii4_3',ii]+
         state['inc_year_gii_3',ii])/(state['pop_by4age3',ii])) ,
  
  ( (  state['inc_year_gi3_4',ii]+
         state['inc_year_gi_4',ii]+
         state['inc_year_gii4_4',ii]+
         state['inc_year_gii_4',ii])/(state['pop_by4age4',ii])) ,
  
  ( (  state['inc_year_gi3_5',ii]+
         state['inc_year_gi_5',ii]+
         state['inc_year_gii4_5',ii]+
         state['inc_year_gii_5',ii])/(state['pop_by4age5',ii])) 
)*1000



iid2_all<-iid2_plot_func(modelled_irate, data_plots$data_iid2.c4,   "all"," ", "F",nice_cols[8])





# GII4 prevalence in children  --------------------------------------------


x_d <- c(1, 2, 3, 4, 5, 6)-0.28 

#ii<-which(!is.na(data_all$sero1))
ii<-data_all$time[which(!is.na(data_all$sero1))]-1 # model strats from 0 and to match by week output

observed_size<-c(
  103*10,
  107*10,
  121*10,
  124*10,
  122*10,
  109*10
)
sero_obs<-(c(data_all$sero1[ii], data_all$sero2[ii], data_all$sero3[ii], data_all$sero4[ii],
             data_all$sero5[ii], data_all$sero6[ii])/observed_size)*1e2

df_d <- data.frame(
  x<-x_d,
  sero= data_plots$dfsero$mean*100,
  low=data_plots$dfsero$low*100,
  up=data_plots$dfsero$up*100)

fac<-c(2,3,4,5,6,7)


#id<-which(tt%in%which(!is.na(data$sero1)))

sero_model<-rbind(
  state['seroprev1.2',ii] ,
  state['seroprev2.3',ii] ,
  state['seroprev3.4',ii] ,
  state['seroprev4.5',ii] ,
  state['seroprev5.6',ii] ,
  state['seroprev6.7',ii])*100



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
       y = "GII.4\n prevalence (%)") +
  theme_classic() +
  
  ylim(0, max(df_d$up)*1.2) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10),
    plot.tag = element_text(size = 10),
    plot.tag.position = c(0.15, 0.91),
    panel.border = element_rect(colour = "black", fill=NA)
  )




# SGSS age reporting fractions --------------------------------------------


#ii<-which(!is.na(data_all$a0_event))
ii<-data_all$time[which(!is.na(data_all$a0_event))]-1 # model strats from 0 and to match by week output



a0_model<-((state['cumm_incday_gi3_1',ii]+
              state['cumm_incday_gi_1',ii]+
              state['cumm_incday_gii4_1',ii]+
              state['cumm_incday_gii_1',ii]) /
             (    
               state['cumm_incday_gii4_1',ii]+
                 state['cumm_incday_gii4_2',ii]+
                 state['cumm_incday_gii4_3',ii]+
                 state['cumm_incday_gii4_4',ii]+
                 state['cumm_incday_gii_1',ii]+
                 state['cumm_incday_gii_2',ii]+
                 state['cumm_incday_gii_3',ii]+
                 state['cumm_incday_gii_4',ii]+
                 state['cumm_incday_gi3_1',ii]+
                 state['cumm_incday_gi3_2',ii]+
                 state['cumm_incday_gi3_3',ii]+
                 state['cumm_incday_gi3_4',ii]+
                 state['cumm_incday_gi_1',ii]+
                 state['cumm_incday_gi_2',ii]+
                 state['cumm_incday_gi_3',ii]+
                 state['cumm_incday_gi_4',ii]))


a5_model<-((state['cumm_incday_gi3_2',ii]+
              state['cumm_incday_gi_2',ii]+
              state['cumm_incday_gii4_2',ii]+
              state['cumm_incday_gii_2',ii]) /
             (    
               state['cumm_incday_gii4_1',ii]+
                 state['cumm_incday_gii4_2',ii]+
                 state['cumm_incday_gii4_3',ii]+
                 state['cumm_incday_gii4_4',ii]+
                 state['cumm_incday_gii_1',ii]+
                 state['cumm_incday_gii_2',ii]+
                 state['cumm_incday_gii_3',ii]+
                 state['cumm_incday_gii_4',ii]+
                 state['cumm_incday_gi3_1',ii]+
                 state['cumm_incday_gi3_2',ii]+
                 state['cumm_incday_gi3_3',ii]+
                 state['cumm_incday_gi3_4',ii]+
                 state['cumm_incday_gi_1',ii]+
                 state['cumm_incday_gi_2',ii]+
                 state['cumm_incday_gi_3',ii]+
                 state['cumm_incday_gi_4',ii]))

a15_model<-((state['cumm_incday_gi3_3',ii]+
               state['cumm_incday_gi_3',ii]+
               state['cumm_incday_gii4_3',ii]+
               state['cumm_incday_gii_3',ii]) /
              (    
                state['cumm_incday_gii4_1',ii]+
                  state['cumm_incday_gii4_2',ii]+
                  state['cumm_incday_gii4_3',ii]+
                  state['cumm_incday_gii4_4',ii]+
                  state['cumm_incday_gii_1',ii]+
                  state['cumm_incday_gii_2',ii]+
                  state['cumm_incday_gii_3',ii]+
                  state['cumm_incday_gii_4',ii]+
                  state['cumm_incday_gi3_1',ii]+
                  state['cumm_incday_gi3_2',ii]+
                  state['cumm_incday_gi3_3',ii]+
                  state['cumm_incday_gi3_4',ii]+
                  state['cumm_incday_gi_1',ii]+
                  state['cumm_incday_gi_2',ii]+
                  state['cumm_incday_gi_3',ii]+
                  state['cumm_incday_gi_4',ii]))

a65_model<-((  state['cumm_incday_gi3_4',ii]+
                 state['cumm_incday_gi_4',ii]+
                 state['cumm_incday_gii4_4',ii]+
                 state['cumm_incday_gii_4',ii]) /
              (    
                state['cumm_incday_gii4_1',ii]+
                  state['cumm_incday_gii4_2',ii]+
                  state['cumm_incday_gii4_3',ii]+
                  state['cumm_incday_gii4_4',ii]+
                  state['cumm_incday_gii_1',ii]+
                  state['cumm_incday_gii_2',ii]+
                  state['cumm_incday_gii_3',ii]+
                  state['cumm_incday_gii_4',ii]+
                  state['cumm_incday_gi3_1',ii]+
                  state['cumm_incday_gi3_2',ii]+
                  state['cumm_incday_gi3_3',ii]+
                  state['cumm_incday_gi3_4',ii]+
                  state['cumm_incday_gi_1',ii]+
                  state['cumm_incday_gi_2',ii]+
                  state['cumm_incday_gi_3',ii]+
                  state['cumm_incday_gi_4',ii]))






age_model<-cbind(a0_model,a5_model, a15_model,a65_model )*100


x_d <- c(1,2,3,4)-0.28 # bin x axis positions

df_d <- data.frame(
  x = x_d,
  age = data_plots$agg_age[,2]*100,
  low = data_plots$agg_age[,1]*100,
  up = data_plots$agg_age[,3]*100
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
  geom_violin(
    data = df_m,
    aes(x = variable, y = value, fill = "Posterior Density"),
    alpha = .4, 
    trim = FALSE
  ) +
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
    axis.title = element_text(size = 10),
    plot.tag = element_text(size = 10),
    plot.tag.position = c(0.2, 0.91),
    panel.border = element_rect(colour = "black", fill=NA)
  )





# Get plots together  -----------------------------------------------------
library(grid)
blank<-grid.rect(gp=gpar(col="white"))


#windows()
gridExtra::grid.arrange(
  fits_sgss,
  fits_sgss_avg,
  blank,
  fits_bystrain,
  iid2_all,
  fits_age,
  fits_sero,
  layout_matrix = rbind(c(1 ,1),
                        c(2 ,3),
                        c(4, 5),
                        c(6, 7)))





