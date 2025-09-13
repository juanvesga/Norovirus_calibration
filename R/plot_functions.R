# -------------------------------------------------------------------------
# Plotting functions-------------------------------------------------------
# -------------------------------------------------------------------------

## trend of reported by strain
plot_reported_strain<-function(strain,dat,num){
  
  df_d<-dat
 
  browser()
  df_s <- as.data.frame(
    rowQuantiles((strain),
                 probs = c(0.025, 0.5, 0.975))
  )
  df_s$date<-df_d$x
  
  
  df_y<-as.data.frame((strain))
  column_names<-names(df_y)
  df_y$date<-df_d$x
  
  
  df_m <- gather(df_y, itr, value, V1:column_names[length(column_names)])
  #data_long$strain<-"GI.3"
  
  
  cols<-c("dodgerblue", "#cc0044","#FFB400","#00A04B","grey")
  
  
  start_date <- df_d$x[1]
  end_date <- df_d$x[length(df_d$x)] 
  
  tag_lab<-c("B","C","D","E")
  
  ylabs<-c("Weekly \n cases reported"," ","Weekly \n cases reported"," ")
  
  fits_sgss <- ggplot() +
    geom_line(data = df_m, aes(x = date, y=value, group=itr),
              col = cols[5] , alpha = 0.2, lwd = 0.01) +
    geom_ribbon(data=df_s,aes(x=date,ymin = `2.5%`,ymax = `97.5%`), fill=cols[5], alpha=0.2)+
    geom_line(data = df_s, aes(x = date, y=`50%`), col = "black",lwd = 0.8) +
    geom_bar(data = df_d, aes(x = x, y = y), stat="identity", #width = 0.9, 
             fill=cols[num], alpha=0.5) +
    labs(tag =tag_lab[num], x = "", y = ylabs[num]) +
    theme_classic() +
    ylim(c(0,max(df_s$`97.5%`)*1.1))+
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
      plot.tag.position = c(0.25, 0.95),
      panel.border = element_rect(colour = "black", fill=NA)
    )
  
  return(fits_sgss)
  
}

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
    geom_violinhalf(
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
         y =paste(title_string , "Incidence \n per 1k person-year",sep=" ")) +
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
  
  return(iid2_plot)
  
}
