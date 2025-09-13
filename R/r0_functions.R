get_R0<-function(pars,params,m,time_vec,get_addresses,strain){
  
  
  
  #fixed model parameters   
  rr_inf_asymp <- 0.05 
  epsilon <- 1.0   # incubation (E->I)
  theta_5 <- 2.5   # duration symptoms in under 5 (I->A)
  theta_5p <- 1.5   # duration symptoms in over 5 
  sigma <- 15 #(A->R)
  mu<-params$mu
  #w2  <- -2
  nonsecretor<-0
  C <- m
  
  
  # CM<-list()
  # CM$cmx1<-params$cmx_1
  # CM$cmx2<-params$cmx_2
  # CM$cmx3<-params$cmx_3
  # CM$cmx4<-params$cmx_4
  # CM$cmx5<-params$cmx_5
  # CM$cmx6<-params$cmx_6
  # CM$cmx7<-params$cmx_7
  # CM$cmx8<-params$cmx_8
  # CM$cmx9<-params$cmx_9
  # 
  # CM$c<-params$contact$matrix
  # CM$c_holi<-params$contact_holi$matrix
  
  
  
  aduRR     <-unlist(pars["aduRR"])
  beta_gi3  <-unlist(pars["beta_1"])    
  beta_gi   <-unlist(pars["beta_2"]) 
  beta_gii4 <-unlist(pars["beta_3"] )  
  beta_gii  <-unlist(pars["beta_4"])
  crossp_GI <-unlist(pars[["crossp_GI"]])
  crossp_GII<-unlist(pars[["crossp_GII"]])
  #w1        <-unlist(pars["w1_1"])
  
  # Sequence of new strain 1 infections after acquiring otehr strains (i.e, 2,3,4)
  gps<-list(
    states=c("E","I","A"),
    age   =c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
  
  state_infectivity<-c(rr_inf_asymp,1,rr_inf_asymp) # compartment relative infectivity
  age_infectivity<-c(1:nrow(C))*0 + 1 # age relative infectivity
  age_infectivity[5:nrow(C)]<-aduRR
  
  theta<-c(1:nrow(C))*0 + theta_5
  theta[6:nrow(C)]<-theta_5p
  
  
  # Select strain-specific parameters 
  if(strain=="GI3"){
    q= beta_gi3
    cross_immm = crossp_GI
  }else if(strain=="GI"){
    q= beta_gi
    cross_immm = crossp_GI
  }else if(strain=="GII4"){
    q= beta_gii4
    nonsecretor<-0.2
    cross_immm = crossp_GII
  }else if (strain=="GII"){
    q= beta_gii
    cross_immm = crossp_GII
  }
  
  
  # Model states
  
  
  groups<-list(gps$states,gps$age)
  i<-list()
  s<-list()
  ref <- get_addresses(groups, i, s, 0)
  
  nage<-length(gps$age)
  nstates<-length(gps$states)
  
  Transition<-matrix(0,nrow = nstates*nage,ncol = nstates*nage)
  
  
  #Transitions
  
  for (ia in 1:length(gps$age)){
    age_i <- gps$age[ia]
    
    gi <- function(sta) ref$i[[sta]][[age_i]]
    E   <- gi('E')
    I   <- gi('I')
    A   <- gi('A')
    
    
    # Transitions in X1
    source <- E 
    destin <- E 
    rate <- 1/epsilon + mu[ia] 
    Transition[destin,source] <-Transition[destin,source] + rate
    
    #  I
    source <- I 
    destin <- I 
    rate <- 1/theta[ia]+ mu[ia]
    Transition[destin,source] <-Transition[destin,source] + rate
    
    #  A
    source <- A 
    destin <- A 
    rate <- 1/sigma+ mu[ia]
    Transition[destin,source] <-Transition[destin,source] + rate
    
    # E -> I 
    source <- E 
    destin <- I 
    rate <- -1/epsilon
    Transition[destin,source] <-Transition[destin,source] + rate
    
    
    # I -> A 
    source <- I 
    destin <- A 
    rate <- - 1/theta[ia]
    Transition[destin,source] <-Transition[destin,source] + rate
    
  }
  
  
  
  Transmission<- matrix(0,nrow = nstates*nage,ncol = nstates*nage)
  Transmission_cross<- matrix(0,nrow = nstates*nage,ncol = nstates*nage)
  Transmission_re<- matrix(0,nrow = nstates*nage,ncol = nstates*nage)
  
  
  for (ia in 1:length(gps$age)){
    age_i <- gps$age[ia]
    
    gi <- function(state_i) ref$i[[state_i]][[age_i]]
    E_i   <- gi('E')
    I_i   <- gi('I')
    A_i   <- gi('A')
    
    for (ja in 1:length(gps$age)){
      age_j <- gps$age[ja]
      
      gj <- function(state_j) ref$i[[state_j]][[age_j]]
      E_j   <- gj('E')
      I_j   <- gj('I')
      A_j   <- gj('A')
      
      ## Transmission fully susceptible 
      
      
      source <- E_j
      infected <- E_i 
      rate <- q*C[ia,ja]*state_infectivity[1] * age_infectivity[ja]  
      Transmission[infected,source] <- Transmission[infected,source] + rate
      
      source <- I_j
      infected <- E_i 
      rate <- q*C[ia,ja]*state_infectivity[2] * age_infectivity[ja] 
      Transmission[infected,source] <- Transmission[infected,source] + rate
      
      source <- A_j
      infected <- E_i 
      rate <- q*C[ia,ja]*state_infectivity[3] * age_infectivity[ja] 
      Transmission[infected,source] <- Transmission[infected,source] + rate
      
      
      ## Transmission with cross protection
      source <- E_j
      infected <- E_i 
      rate <- q*C[ia,ja]*state_infectivity[1] * age_infectivity[ja]  * (1-cross_immm)
      Transmission_cross[infected,source] <- Transmission_cross[infected,source] + rate
      
      source <- I_j
      infected <- E_i 
      rate <- q*C[ia,ja]*state_infectivity[2] * age_infectivity[ja] * (1-cross_immm)
      Transmission_cross[infected,source] <- Transmission_cross[infected,source] + rate
      
      source <- A_j
      infected <- E_i 
      rate <- q*C[ia,ja]*state_infectivity[3] * age_infectivity[ja]  * (1-cross_immm)
      Transmission_cross[infected,source] <- Transmission_cross[infected,source] + rate
      
      source <- E_j
      infected <- A_i 
      rate <- q*C[ia,ja]*state_infectivity[1] * age_infectivity[ja]  * (cross_immm)
      Transmission_cross[infected,source] <- Transmission_cross[infected,source] + rate
      
      source <- I_j
      infected <- A_i 
      rate <- q*C[ia,ja]*state_infectivity[2] * age_infectivity[ja] * (cross_immm)
      Transmission_cross[infected,source] <- Transmission_cross[infected,source] + rate
      
      source <- A_j
      infected <- A_i 
      rate <- q*C[ia,ja]*state_infectivity[3] * age_infectivity[ja]  * (cross_immm)
      Transmission_cross[infected,source] <- Transmission_cross[infected,source] + rate
      
      ############### asymptomatic reinfection
      source <- E_j
      infected <- A_i 
      rate <- q*C[ia,ja]*state_infectivity[1] * age_infectivity[ja] 
      Transmission_re[infected,source] <- Transmission_re[infected,source] + rate
      
      source <- I_j
      infected <- A_i 
      rate <- q*C[ia,ja]*state_infectivity[2] * age_infectivity[ja] 
      Transmission_re[infected,source] <- Transmission_re[infected,source] + rate
      
      source <- A_j
      infected <- A_i 
      rate <- q*C[ia,ja]*state_infectivity[3] * age_infectivity[ja] 
      Transmission_re[infected,source] <- Transmission_re[infected,source] + rate
      
      
    }
  }
  
  
  
  
  
  library(matlib)
  
  V_inv<-solve(Transition)
  
  ########### Full susceptibility infection
  KL<-Transmission%*%V_inv
  nz<-which(rowSums(KL)!=0)
  E<-matrix(0,nage*nstates, ncol=length(nz))
  
  for (ii in 1:length(nz)){
    E[nz[ii],ii]<-1
  }
  
  # Small K
  K <- t(E)%*%KL%*%E
  NGM<-K
  R0_symp=max(Re(eigen(NGM)$values))
  
  ########### cross-protection infection
  KL<-Transmission_cross%*%V_inv
  nz<-which(rowSums(KL)!=0)
  E<-matrix(0,nage*nstates, ncol=length(nz))
  
  for (ii in 1:length(nz)){
    E[nz[ii],ii]<-1
  }
  
  # Small K
  K <- t(E)%*%KL%*%E
  NGM<-K
  R0_cross=max(Re(eigen(NGM)$values))
  
  
  
  ########### reinfection R
  KL<-Transmission_re%*%V_inv
  nz<-which(rowSums(KL)!=0)
  E<-matrix(0,nage*nstates, ncol=length(nz))
  for (ii in 1:length(nz)){
    E[nz[ii],ii]<-1
  }
  
  # Small K
  K <- t(E)%*%KL%*%E
  NGM<-K
  R0_re=max(Re(eigen(NGM)$values))
  
  tt<-seq(1,time_vec,1)
  
  seasonality <- (1 + 0.15* cos(2 * pi * (tt - 10) /365))
  
  
  R0_full <- R0_re + R0_symp + R0_cross
  Rs_re <- (R0_re)*seasonality 
  Rs_symp <- (R0_symp)*seasonality 
  Rs_re <- (R0_re)*seasonality
  Rs_cross<- (R0_cross)*seasonality 
  season <-  seasonality
  
  
  
  res<-list()
  res$R0_symp <- mean(R0_symp)
  res$R0_cross <- mean(R0_cross)
  res$R0_re <- mean(R0_re)
  res$R0_full <- mean(R0_re + R0_symp + R0_cross)
  res$Rs_re <- mean(R0_re)*seasonality 
  res$Rs_symp <- mean(R0_symp)*seasonality 
  res$Rs_re <- mean(R0_re)*seasonality 
  res$Rs_cross<- mean(R0_cross)*seasonality 
  res$season <-  season
  
  return(res)
}





get_addresses<- function (groups, i, s, lim){
  
  
  # Initiate any sets not so far covered by s
  if  (!length(s)){
    fnames = {}
  } else {
    fnames <- names(s)
  }
  
  for (ig in  1:length(groups)){
    gp <- groups[[ig]]
    for (ig2 in 1:length(gp)){
      if ( !(gp[ig2] %in% names(fnames))){  
        s[[gp[ig2]]]<-{}
      }
    }
  }
  
  
  if (length(groups) == 1){
    gp1 <- groups[[1]]
    for (ig1 in 1:length(gp1)){
      lim <- lim+1
      i[[gp1[ig1]]] <- lim
      s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
      
    }
  }
  
  
  if (length(groups) == 2){
    gp1 = groups[[1]] 
    gp2 = groups[[2]]
    for (ig1 in 1:length(gp1)){
      i[[gp1[ig1]]]<-list()
      for (ig2 in 1:length(gp2)){
        lim <- lim+1
        i[[gp1[ig1]]][[gp2[ig2]]]<-lim
        s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
        s[[gp2[ig2]]] <- c(s[[gp2[ig2]]], lim)
        
      }
    }
  }
  
  if (length(groups) == 3){
    gp1 = groups[[1]] 
    gp2 = groups[[2]]
    gp3 = groups[[3]]
    for (ig1 in 1:length(gp1)){
      i[[gp1[ig1]]]<-list()
      for (ig2 in 1:length(gp2)){
        i[[gp1[ig1]]][[gp2[ig2]]]<-list()
        for (ig3 in 1:length(gp3)){
          lim <- lim+1
          i[[gp1[ig1]]][[gp2[ig2]]][[gp3[ig3]]]<-lim
          
          s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
          s[[gp2[ig2]]] <- c(s[[gp2[ig2]]], lim)
          s[[gp3[ig3]]] <- c(s[[gp3[ig3]]], lim)
          
        }
      }
    }
  }
  
  if (length(groups) == 4){
    gp1 = groups[[1]] 
    gp2 = groups[[2]]
    gp3 = groups[[3]]
    gp4 = groups[[4]]
    for (ig1 in 1:length(gp1)){
      i[[gp1[ig1]]]<-list()
      for (ig2 in 1:length(gp2)){
        i[[gp1[ig1]]][[gp2[ig2]]]<-list()
        for (ig3 in 1:length(gp3)){
          i[[gp1[ig1]]][[gp2[ig2]]][[gp3[ig3]]]<-list()
          for (ig4 in 1:length(gp4)){
            lim <- lim+1
            i[[gp1[ig1]]][[gp2[ig2]]][[gp3[ig3]]][[gp4[ig4]]]<-lim
            
            s[[gp1[ig1]]] <- c(s[[gp1[ig1]]], lim)
            s[[gp2[ig2]]] <- c(s[[gp2[ig2]]], lim)
            s[[gp3[ig3]]] <- c(s[[gp3[ig3]]], lim)
            s[[gp4[ig4]]] <- c(s[[gp4[ig4]]], lim)
            
          }
        }
      }
    }
  }
  
  i$nstates<-lim
  return(list(s=s, i=i))
}
