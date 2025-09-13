# index 1 is GI3 infection
# index 2 is other Gi infection
# index 3 is GII4 infection
# index 4 is other GII infection
# For states E, I and A the first number is the active infecyive strain
# the following numbers are the carrying immunity
# R compartments show carrying immunity in ascending order

## Equations for transitions between compartments by age group
###############################################################
# Not infected
noro_model <- odin({
 
  
  update(G[]) <-  
    G[i] + 
    n_bG[i] - 
    n_ageoG[i] + 
    n_ageiG[i] - 
    n_muG[i]  
  
  update(M[]) <-  
    M[i] + 
    n_bM[i] - 
    n_ageoM[i] + 
    n_ageiM[i] - 
    n_MS[i] - 
    n_muM[i]  
  
  update(S[]) <-
    S[i] + 
    n_bS[i] - 
    n_ageoS[i] + 
    n_ageiS[i] + 
    n_MS[i] + 
    sum(n_R_wane[i,])-
    sum(n_S_E[i,]) -
    n_muS[i] 
  
  
  dim(n_S_E)<-c(N_age,N_strain)
  # Infection # 1: age x 4 dimensions
  update(E[,]) <-
    E[i,j]- 
    n_ageoE[i,j] + 
    n_ageiE[i,j] + 
    n_R_E[i,j] +
    n_S_E[i,j]  - 
    n_E_I[i,j] - 
    n_muE[i,j] 
  
  dim(n_E_I)<-c(N_age,N_strain)
  dim(n_R_E)<-c(N_age,N_strain)
  
  update(I[,]) <-   
    I[i,j] - 
    n_ageoI[i,j] + 
    n_ageiI[i,j] + 
    n_E_I[i,j] -
    n_I_A[i,j] - 
    n_muI[i,j] 
  
  dim(n_I_A)<-c(N_age,N_strain)
  
  update(A[,]) <-  
    A[i,j]- 
    n_ageoA[i,j] + 
    n_ageiA[i,j] + 
    n_I_A[i,j] + 
    n_R_A[i,j] - 
    n_A_R[i,j] - 
    n_muA[i,j]
  
  dim(n_R_A)<-c(N_age,N_strain)
  dim(n_A_R)<-c(N_age,N_strain)
  
  update(R[,]) <- 
    R[i,j] - 
    n_ageoR[i,j] + 
    n_ageiR[i,j] + 
    n_A_R[i,j]   - 
    n_R_wane[i,j]-
    sum(n_R_tot[i,j,]) -
    # n_R_E[i,j] -
    # n_R_A[i,j] -
    n_muR[i,j]
  
  dim(n_R_wane)<-c(N_age,N_strain)
  
 
  
  
  ####### Outputs
  #########################
  reported_fac[1]<-1/repfac_0
  reported_fac[2]<-1/repfac_5
  reported_fac[3]<-1/repfac_15
  reported_fac[4]<-1/repfac_65p
  dim(reported_fac)<-4
  
  reported_fac_long[1:4]<-1/repfac_0
  reported_fac_long[5:8]<-1/repfac_5
  reported_fac_long[9:13]<-1/repfac_15
  reported_fac_long[N_age]<-1/repfac_65p
  dim(reported_fac_long)<-N_age
  
  reported_geno[1]<-geno_frac*1/repfac_0
  reported_geno[2]<-geno_frac*1/repfac_5
  reported_geno[3]<-geno_frac*1/repfac_15
  reported_geno[4]<-geno_frac*1/repfac_65p
  dim(reported_geno)<-4
  
  
  # Daily infections incidence
  update(infections_day_gi3[]) <- 
    if (i==1)
      (infections_day_gi3[i]+
         sum(n_E_I[1:4,1]) * reported_fac[i] 
      ) else
        if(i==2) 
          (infections_day_gi3[i]+
             sum(n_E_I[5:8,1]) * reported_fac[i]  
          ) else
            if(i==3) 
              (infections_day_gi3[i]+
                 sum(n_E_I[9:13,1]) * reported_fac[i]  
              ) else
                (infections_day_gi3[i]+
                   sum(n_E_I[14,1]) * reported_fac[i] )
  
  
  update(infections_day_gi[]) <- 
    if (i==1)
      (infections_day_gi[i]+
         sum(n_E_I[1:4,2]) * reported_fac[i]   
      ) else
        if(i==2) 
          (infections_day_gi[i]+
             sum(n_E_I[5:8,2]) * reported_fac[i]  
          ) else
            if(i==3) 
              (infections_day_gi[i]+
                 sum(n_E_I[9:13,2]) * reported_fac[i] 
              ) else
                (infections_day_gi[i]+
                   sum(n_E_I[14,2]) * reported_fac[i]  )
  
  
  update(infections_day_gii4[]) <- 
    if (i==1)
      (infections_day_gii4[i]+
         sum(n_E_I[1:4,3]) * reported_fac[i]  
      ) else
        if(i==2) 
          (infections_day_gii4[i]+
             sum(n_E_I[5:8,3]) * reported_fac[i]  
          ) else
            if(i==3) 
              (infections_day_gii4[i]+
                 sum(n_E_I[9:13,3]) * reported_fac[i]   
              ) else
                (infections_day_gii4[i]+
                   sum(n_E_I[14,3]) * reported_fac[i] )
  
  
  
  update(infections_day_gii[]) <- 
    if (i==1)
      (infections_day_gii[i]+
         sum(n_E_I[1:4,4]) * reported_fac[i]  
      ) else
        if(i==2) 
          (infections_day_gii[i]+
             sum(n_E_I[5:8,4]) * reported_fac[i]    
          ) else
            if(i==3) 
              (infections_day_gii[i]+
                 sum(n_E_I[9:13,4]) * reported_fac[i]    
              ) else
                (infections_day_gii[i]+
                   sum(n_E_I[14,4]) * reported_fac[i]  )
  
  
  # Weekly reported cases (match sgss)
  
  update(reported_wk[]) <- reported_wk[i] + 
    (sum(n_E_I[i,]) * reported_fac_long[i] )  
  
  # update(death_reported_wk_all[]) <- death_reported_wk_all[i] +
  #   (reported_wk[i] * p_cfr[i])
  # 
  # update(death_reported_wk) <- death_reported_wk +
  #   sum(death_reported_wk_all)
  
  update(reported_wk_gi3) <- reported_wk_gi3 + 
    (sum(n_E_I[1:4,1])* reported_geno[1] ++ 
       sum(n_E_I[5:8,1])* reported_geno[2] +
       sum(n_E_I[9:13,1])* reported_geno[3] +
       sum(n_E_I[14,1])* reported_geno[4] ) 
  
  
  
  update(reported_wk_gi) <- reported_wk_gi + 
    (sum(n_E_I[1:4,2] )* reported_geno[1]  +
       sum(n_E_I[5:8,2] )* reported_geno[2] +
       sum(n_E_I[9:13,2] )* reported_geno[3] +
       sum(n_E_I[14,2] )* reported_geno[4]  ) 
  
  
  
  update(reported_wk_gii4) <- reported_wk_gii4 + 
    (sum(n_E_I[1:4,3])* reported_geno[1] +
       sum(n_E_I[5:8,3])* reported_geno[2] +
       sum(n_E_I[9:13,3])* reported_geno[3] +
       sum(n_E_I[14,3])* reported_geno[4] ) 
  
  
  
  update(reported_wk_gii) <- reported_wk_gii + 
    (sum(n_E_I[1:4,4])* reported_geno[1] +
       sum(n_E_I[5:8,4])* reported_geno[2] +
       sum(n_E_I[9:13,4])* reported_geno[3] +
       sum(n_E_I[14,4])* reported_geno[4] ) 
  
  
  
  ## Incidence by Year and strain
  ###############
  
  
  update(inc_year_gi3[]) <- 
    if(i==1) 
      inc_year_gi3[i]+
    (sum(n_E_I[1,1])) else if(i==2)  
         inc_year_gi3[i]+
    sum(n_E_I[2:4,1])  else if(i==3)  
      inc_year_gi3[i]+
    sum(n_E_I[5:8,1])  else if (i==4) 
      inc_year_gi3[i]+
    sum(n_E_I[9:13,1])   else
      inc_year_gi3[i]+
    sum(n_E_I[14,1])
  
  
  
  
  
  update(inc_year_gi[]) <-  
    if(i==1) 
      inc_year_gi[i]+
    sum(n_E_I[1,2] ) 
     else if(i==2)  
      inc_year_gi[i]+
    sum(n_E_I[2:4,2] ) 
    else if(i==3) 
      inc_year_gi[i]+
    sum(n_E_I[5:8,2] ) 
    else if (i==4) 
      inc_year_gi[i]+
    sum(n_E_I[9:13,2] )  else
      inc_year_gi[i]+
    sum(n_E_I[14,2] ) 
  
  
  update(inc_year_gii4[]) <- 
    if(i==1) 
      inc_year_gii4[i]+
    sum(n_E_I[1,3])   else if(i==2)  
      inc_year_gii4[i]+
    sum(n_E_I[2:4,3])  else if(i==3)  
      inc_year_gii4[i]+
    sum(n_E_I[5:8,3])  else if (i==4) 
      inc_year_gii4[i]+
    sum(n_E_I[9:13,3])   else
      inc_year_gii4[i]+
    sum(n_E_I[14,3])  
  
  update(inc_year_gii[]) <-  
    if(i==1) 
      inc_year_gii[i]+
    sum(n_E_I[1,4])   else if(i==2)  
      inc_year_gii[i]+
    sum(n_E_I[2:4,4])  else if(i==3)  
      inc_year_gii[i]+
    sum(n_E_I[5:8,4])  else if (i==4) 
      inc_year_gii[i]+
    sum(n_E_I[9:13,4])    else
      inc_year_gii[i]+
    sum(n_E_I[14,4])  
  
  
  ## compute at risk population ( match IDD2 data)
  update(pop_by4age[]) <-
    (if(i==1)  sum(N_byage[1]) else
      if(i==2)  sum(N_byage[2:4])  else 
        if(i==3)  sum(N_byage[5:8])  else 
          if (i==4)  sum(N_byage[9:13])  else
            sum(N_byage[14]) ) 
  
  ## compute prevalence og GII4
  
  update(seroprev_gii4[]) <- prev_byage_gii4[i]/N_byage[i]
  
  
  
  # Individual probabilities of transition: -------------------------------------------------------
  
  p_mu[]     <- 1 - exp(-mu[i] * dt) # mortality
  p_aging[]  <- 1 - exp(-aging_vec[i] * dt)
  p_MS[]     <- 1 - exp(-(1/maternalAB) * dt)  # M to S
  p_SE[,]    <- 1 - exp(- lambda[i,j] * dt) # S to E all
  p_EI       <- 1 - exp(-(1/epsilon) * dt) # E to I
  p_IA_5     <- 1 - exp(-(1/theta_5) * dt) # I to A
  p_IA_5p    <- 1 - exp(-(1/theta_5p) * dt) # I to A
  p_AR       <- 1 - exp(-(1/sigma) * dt) # A to R
  p_wane       <- 1 - exp(- (1/(imm_yr*365)) * dt) # R to S
  
  N_byage[]<-(
      G[i] +
      M[i] +  
      S[i] + 
      sum(E[i,])+
      sum(I[i,])+
      sum(A[i,])+
      sum(R[i,]))
  
  
  
  prev_byage_gii4[]<-(M[i] +
                        I[i,3]+
                        A[i,3]+
                        R[i,3])
  
  prev_byage_all[]<-(
    M[i]+
      sum(I[i,])+
      sum(A[i,])+
      sum(R[i,]))
  
  
  
  
  dim(N_byage)<-N_age
  dim(prev_byage_gii4)<-N_age
  dim(prev_byage_all)<-N_age
  
  N <- sum(N_byage)
  
  prev   <- sum(prev_byage_all)/N  
  
  
  
  
  ## Force of infection
  ###########################
  
  
  school_switcher <- interpolate(school_time,school_value,"constant")
  
  dim(school_time, school_value) <- parameter(rank = 1)
  
  comix_switcher <- interpolate(comix_time,comix_value,"constant")
  
  dim(comix_time, comix_value) <- parameter(rank = 1)
  
  
  
  contact_matrix[,]<- if (comix_switcher >0)( 
    if (comix_switcher==1)
      cmx_1[i,j] else 
        if (comix_switcher==2)
          cmx_2[i,j] else 
            if (comix_switcher==3)
              cmx_3[i,j] else 
                if (comix_switcher==4)
                  cmx_4[i,j] else 
                    if (comix_switcher==5)
                      cmx_5[i,j] else 
                        if (comix_switcher==6)
                          cmx_6[i,j] else 
                            if (comix_switcher==7)
                              cmx_7[i,j] else 
                                if (comix_switcher==8)
                                  cmx_8[i,j] else 
                                    cmx_9[i,j] 
  ) else (
    school_switcher*m[i, j]  + (1-school_switcher)*m_holi[i, j])
  
  
  
  
  # Age infectivity
  dim(age_RR)<-N_age
  
  age_RR[]<-if(i<8) 1 else aduRR
  
  
  # infectives GI3
  c1_ij[] <- (I[i,1] + 
                rr_inf_asymp * (A[i,1] +
                                  E[i,1]) ) * age_RR[i] 
  
  
  
  # infectives GI
  c2_ij[] <- (I[i,2]+
                rr_inf_asymp * (A[i,2]+
                                  E[i,2]))   * age_RR[i] 
  
  
  # infectives GII4
  c3_ij[] <- (I[i,3]+
                rr_inf_asymp *(A[i,3]+
                                 E[i,3]) )  * age_RR[i] 
  
  
  
  
  c4_ij[] <- (I[i,4]+
                rr_inf_asymp *(A[i,4]+
                                 E[i,4]) )   * age_RR[i] 
  
  dim(c1_ij) <-  N_age
  dim(c2_ij) <-  N_age
  dim(c3_ij) <-  N_age
  dim(c4_ij) <-  N_age
  s_ij_1[,]<- contact_matrix[i,j]  * c1_ij[j]
  s_ij_2[,]<- contact_matrix[i,j]  * c2_ij[j]
  s_ij_3[,]<- contact_matrix[i,j]  * c3_ij[j]
  s_ij_4[,]<- contact_matrix[i,j]  * c4_ij[j]
  
  
  dim(s_ij_1) <- c(N_age, N_age)
  dim(s_ij_2) <- c(N_age, N_age)
  dim(s_ij_3) <- c(N_age, N_age)
  dim(s_ij_4) <- c(N_age, N_age)
  
  # Note the -20 is to rescale the parameters and allow calibration pof negative values of lag with a log transformation
  season_func<-(1 + season_amp* cos(2 * pi * (time - (season_lag-20)) /365))
  
  
  beta_1_t <- beta_1 * season_func#(1 + 0.15* cos(2 * pi * (time - 0) /365))
  beta_2_t <- beta_2 * season_func#(1 + 0.15* cos(2 * pi * (time - 0) /365))
  beta_3_t <- beta_3 * season_func#(1 + 0.15* cos(2 * pi * (time - 0) /365))
  beta_4_t <- beta_4 * season_func#(1 + 0.15* cos(2 * pi * (time - 0) /365))
  
  lambda_1[] <-  beta_1_t   * sum(s_ij_1[i, ])
  lambda_2[] <-  beta_2_t   * sum(s_ij_2[i, ])
  lambda_3[] <-  beta_3_t   * sum(s_ij_3[i, ])
  lambda_4[] <-  beta_4_t   * sum(s_ij_4[i, ])
  
  lambda[,1]  <-lambda_1[i] 
  lambda[,2]  <-lambda_2[i] 
  lambda[,3]  <-lambda_3[i] 
  lambda[,4]  <-lambda_4[i] 
  
  
  
  
  ########### Aging numbers ::::::::::::::::::::::::
  ####################################################
  #Age out
  n_ageoG[] <- Binomial(G[i] ,  p_aging[i])
  n_ageoM[] <- Binomial(M[i] ,  p_aging[i])
  n_ageoS[] <- Binomial(S[i] ,  p_aging[i])
  n_ageoE[,] <- Binomial(E[i,j] ,  p_aging[i])
  n_ageoI[,] <- Binomial(I[i,j] ,  p_aging[i])
  n_ageoA[,] <- Binomial(A[i,j] ,  p_aging[i])
  n_ageoR[,] <- Binomial(R[i,j] ,  p_aging[i])
  
  
  
  
  #Age in
  n_ageiG[] <- if (i>1) Binomial(G[i-1] ,p_aging  [i-1]) else 0
  n_ageiM[] <- if (i>1) Binomial(M[i-1] ,p_aging  [i-1]) else 0
  n_ageiS[] <- if (i>1) Binomial(S[i-1] ,p_aging  [i-1]) else 0
  n_ageiE[,] <- if (i>1) Binomial(E[i-1,j] , p_aging[i-1])else 0
  n_ageiI[,] <- if (i>1) Binomial(I[i-1,j] , p_aging[i-1])else 0
  n_ageiA[,] <- if (i>1) Binomial(A[i-1,j] , p_aging[i-1])else 0
  n_ageiR[,] <- if (i>1) Binomial(R[i-1,j] , p_aging[i-1])else 0

  
  ########### Draws from binomial distributions for numbers changing between
  ## compartments:
  ###############################################
  
  ## Binomial draw for mortality
  n_muG[] <- Binomial(G[i] - n_ageoG[i], p_mu[i])
  n_muM[] <- Binomial(M[i] - n_ageoM[i], p_mu[i])
  n_muS[] <- Binomial(S[i] - n_ageoS[i], p_mu[i])
  n_muE[, ] <- Binomial(E[i, j] - n_ageoE[i, j], p_mu[i])
  n_muI[, ] <- Binomial(I[i, j] - n_ageoI[i, j], p_mu[i])
  n_muA[, ] <- Binomial(A[i, j] - n_ageoA[i, j], p_mu[i])
  n_muR[, ] <- Binomial(R[i, j] - n_ageoR[i, j], p_mu[i])
  
  
  
  
  # M transitions -----------------------------------------------------------
  
  
  
  n_MS[]   <- Binomial(M[i] - n_ageoM[i] - n_muM[i], p_MS[i])
  
  
  
  # rel_foi_strain[,]<- if (sum(lambda[i,]) == 0) 0 else(
  #   if (j < N_strain)
  #     max(
  #       1-(1-(min(lambda[i, j] / sum(lambda[i,]),as.numeric(1)))),
  #       as.numeric(0)
  #     ) else
  #       max(min(1-sum(rel_foi_strain[i, 1:3]),as.numeric(1)),as.numeric(0))
  # )
  
  
  ## rel_foi_strain is probability of an infection in age group i, 
  ## being of strain j
  ##
  ## NOTE: the min(x / sum(x), 1) is required here because with floats,
  ## on a GPU, and with fast math, the sum can include sufficient
  ## rounding error that x / sum(x) can be > 1 by a very small amount;
  ## this keeps us bounded correctly.
  
  rel_foi_strain[, ] <-
    (if (sum(lambda[i,]) == 0)
      (if (j == 1) 1 else 0) else
        min(lambda[i, j] / sum(lambda[i,]),
            as.numeric(1)))
  
  
  dim(rel_foi_strain) <- c(N_age, N_strain)
  
  
  # S to E transitions ------------------------------------------------------
  
  
  
  n_S_E_tot[] <- Binomial(S[i] - n_ageoS[i] - n_muS[i], sum(p_SE[i,]))
  
  
  # Multinomial draw (via nested binomial) for multiple strain infection
  n_S_E[, 1] <- Binomial(n_S_E_tot[i] , rel_foi_strain[i, 1])
  n_S_E[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial((n_S_E_tot[i] - sum(n_S_E[i, 1:(j - 1)])) , rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  
  n_S_E[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    n_S_E_tot[i] - sum(n_S_E[i, 1:3]) 
  
  dim(n_S_E_tot) <- N_age
  
  
  ## Transition E to I
  n_E_I[, ] <- Binomial(E[i, j] - n_ageoE[i, j] - n_muE[i, j], p_EI)
  
  
  # Transitions I to A
  n_I_A[, ] <- Binomial(I[i, j] - n_ageoI[i, j] - n_muI[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  
  
  # A to R transitions
  n_A_R[, ] <- Binomial(A[i, j] - n_ageoA[i, j] - n_muA[i, j], p_AR)
  
  
  
  ########## Transitions R 
  cross_map[,] <-if (i==1) ( 
    if(j==2) crossp_GI else 1 ) else if(i==2)(
      if(j==1) crossp_GI else 1) else if(i==3)(
        if(j==4) crossp_GII else 1) else (
          if(j==3) crossp_GII else 1)
  
  dim(cross_map)<-c(N_strain,N_strain)
  
  
  n_R_wane[, ] <- Binomial(R[i, j] - n_ageoR[i, j] - n_muR[i, j], p_wane)
  
  
  
  
  # i = age
  # j = current strain 
  # k = infection with new strain
  
  n_R_tot[, , ] <- Binomial(((R[i, j] - n_ageoR[i, j] - n_muR[i, j] - n_R_wane[i, j] ) ), p_SE[i,k] * cross_map[j,k])
  
  # R[1]
  
  n_R_rel1[, 1] <- Binomial(sum(n_R_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R_tot[i, 1,]) - sum(n_R_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R_tot[i, 1,]) - sum(n_R_rel1[i, 1:3])
  
  n_R_A[,1]  <- n_R_rel1[i,1]
  n_R_A[,2]  <- n_R_rel2[i,2]
  n_R_A[,3]  <- n_R_rel3[i,3]
  n_R_A[,4]  <- n_R_rel4[i,4]
  
  
  # R[2]
  n_R_rel2[, 1] <- Binomial(sum(n_R_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R_tot[i, 2,]) - sum(n_R_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R_tot[i, 2,]) - sum(n_R_rel2[i, 1:3])
  
  
  
  # R[3]
  n_R_rel3[, 1] <- Binomial(sum(n_R_tot[i, 3,]), rel_foi_strain[i, 1])
  n_R_rel3[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R_tot[i, 3,]) - sum(n_R_rel3[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R_rel3[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R_tot[i, 3,]) - sum(n_R_rel3[i, 1:3])
  
  
  # R[4]
  n_R_rel4[, 1] <- Binomial(sum(n_R_tot[i, 4,]), rel_foi_strain[i, 1])
  n_R_rel4[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R_tot[i, 4,]) - sum(n_R_rel4[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R_rel4[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R_tot[i, 4,]) - sum(n_R_rel4[i, 1:3])

  
  
  n_R_E[,1]  <- n_R_rel2[i, 1]  + n_R_rel3[i,1] + n_R_rel4[i,1] 
  n_R_E[,2]  <- n_R_rel1[i, 2]  + n_R_rel3[i,2] + n_R_rel4[i,2]
  n_R_E[,3]  <- n_R_rel4[i, 3]  + n_R_rel1[i,3] + n_R_rel2[i,3]  
  n_R_E[,4]  <- n_R_rel3[i, 4]  + n_R_rel1[i,4] + n_R_rel2[i,4]
  
  
  dim(n_R_tot) <- c(N_age, N_strain,N_strain)
  dim(n_R_rel1) <- c(N_age, N_strain)
  dim(n_R_rel2) <- c(N_age, N_strain)
  dim(n_R_rel3) <- c(N_age, N_strain)
  dim(n_R_rel4) <- c(N_age, N_strain)
  
  
  
  
  ##### Deaths 
  ###################################
  
  n_allDeath<- 
    sum(n_muG)+
    sum(n_muM)+
    sum(n_muS)+
    sum(n_muE)+
    sum(n_muI)+
    sum(n_muA)+
    sum(n_muR)
  
  
  # Births to keep stable population equal to deaths
  n_bG[] <- if (i==1) (n_allDeath*p_nonsecretor) else 0 
  n_bM[] <- if (i==1) (n_allDeath*(1-p_nonsecretor)*prev) else 0 
  n_bS[] <- if (i==1) ((n_allDeath-n_bG[1]-n_bM[1])) else 0 
  
  
  
  ## Initial states:
  initial(G[])<-pop[i]*p_nonsecretor
  initial(M[])<-0
  initial(S[])<-pop[i]*(1-p_nonsecretor)
  initial(E[,])<-0
  initial(I[,])<-if(i==5) 10 else 0
  initial(A[,])<-0
  initial(R[,])<-0
  
  
  # Outputs
  
  initial(inc_year_gii4[], zero_every = 365) <- 0
  initial(inc_year_gii[], zero_every = 365) <- 0
  initial(inc_year_gi3[], zero_every = 365) <- 0
  initial(inc_year_gi[], zero_every = 365) <- 0
  
  initial(pop_by4age[]) <- 0
  initial(seroprev_gii4[]) <-0
  initial(infections_day_gi3[])<-0
  initial(infections_day_gi[])<-0
  initial(infections_day_gii4[])<-0
  initial(infections_day_gii[])<-0
  initial(reported_wk[], zero_every = 7) <- 0
  initial(reported_wk_gi3, zero_every = 7) <- 0
  initial(reported_wk_gi, zero_every = 7) <- 0
  initial(reported_wk_gii4, zero_every = 7) <- 0
  initial(reported_wk_gii, zero_every = 7) <- 0
  
  
  
  
  # Parameters --------------------------------------------------------------
  
  
  # Calibrated Parameters --------------------------------------------------------
  beta_1 <- parameter(0.2)
  beta_2 <- parameter(0.2)
  beta_3 <- parameter(0.2)
  beta_4 <- parameter(0.2)
  maternalAB <- parameter(180)
  aduRR <- parameter(0.1)
  imm_yr <- parameter(5.1)
  repfac_0 <- parameter(287)
  repfac_5 <- parameter(287)
  repfac_15 <- parameter(287)
  repfac_65p <- parameter(287)
  crossp_GI <- parameter(0.05)
  crossp_GII <- parameter(0.05)
  season_lag  <-parameter(10)
  season_amp  <-parameter(0.15)

  
  # Fixed with default--------------------------------------------------------
  epsilon <- parameter(1)
  theta_5 <- parameter(2.5)
  theta_5p <- parameter(1.5)
  sigma <- parameter(15)
  rr_inf_asymp <- parameter(0.05)
  p_nonsecretor <- parameter(0.2)
  geno_frac       <-0.2
  N_strain        <- parameter(4)
  
  
  # Fixed with NO default--------------------------------------------------------
  pop<- parameter()
  mu <- parameter()
  aging_vec <- parameter()
  N_age <- parameter()
  m <- parameter()
  m_holi <- parameter()
  cmx_1 <- parameter()
  cmx_2 <- parameter()
  cmx_3 <- parameter()
  cmx_4 <- parameter()
  cmx_5 <- parameter()
  cmx_6 <- parameter()
  cmx_7 <- parameter()
  cmx_8 <- parameter()
  cmx_9 <- parameter()
  school_time         <- parameter()
  school_value        <- parameter()
  comix_time          <- parameter()
  comix_value         <- parameter()
  
  
  
  # Dimensions--------------------------------------------------------

  dim(aging_vec)<- N_age
  dim(	G	)<-N_age
  dim(	M	)<-N_age
  dim(	S	)<-N_age
  dim(	E	)<-c(N_age,N_strain)
  dim(	I	)<-c(N_age,N_strain)
  dim(	A	)<-c(N_age,N_strain)
  dim(	R	)<-c(N_age,N_strain)


  dim(n_ageiG)<-N_age
  dim(n_ageiM)<-N_age
  dim(n_ageiS)<-N_age
  dim(n_ageiE)<-c(N_age,N_strain)
  dim(n_ageiI)<-c(N_age,N_strain)
  dim(n_ageiA)<-c(N_age,N_strain)
  dim(n_ageiR)<-c(N_age,N_strain)
  
  dim(n_ageoG)<-N_age
  dim(n_ageoM)<-N_age
  dim(n_ageoS)<-N_age
  dim(n_ageoE)<-c(N_age,N_strain)
  dim(n_ageoI)<-c(N_age,N_strain)
  dim(n_ageoA)<-c(N_age,N_strain)
  dim(n_ageoR)<-c(N_age,N_strain)
  
  
  dim(n_muG)<-N_age
  dim(n_muM)<-N_age
  dim(n_muS)<-N_age
  dim(n_muE)<-c(N_age,N_strain)
  dim(n_muI)<-c(N_age,N_strain)
  dim(n_muA)<-c(N_age,N_strain)
  dim(n_muR)<-c(N_age,N_strain)
  
  
  #dim(cumu_inc) <- N_age
  dim(seroprev_gii4) <- 7
  
  
  dim(inc_year_gii4)<- 5
  dim(inc_year_gii)<- 5
  dim(inc_year_gi3)<- 5
  dim(inc_year_gi)<- 5
  dim(pop_by4age)<- 5
  #dim(n_risk) <- N_age
  dim(infections_day_gii4)<- 4
  dim(infections_day_gii)<- 4
  dim(infections_day_gi3)<- 4
  dim(infections_day_gi)<- 4
  dim(reported_wk)<- N_age
  
  
  
  dim(	n_MS	)<-  N_age
  dim(n_bG)   <- N_age
  dim(n_bM)   <- N_age
  dim(n_bS)   <- N_age
  dim(p_mu)   <- N_age
  
  dim(p_aging)<- N_age
  dim(p_MS)   <- N_age
  dim(p_SE)   <- c(N_age,N_strain)
  
  
  dim(pop)<-N_age
  dim(mu) <-N_age
  dim(m) <- c(N_age, N_age)
  dim(m_holi) <- c(N_age, N_age)
  dim(cmx_1) <- c(N_age, N_age)
  dim(cmx_2) <- c(N_age, N_age)
  dim(cmx_3) <- c(N_age, N_age)
  dim(cmx_4) <- c(N_age, N_age)
  dim(cmx_5) <- c(N_age, N_age)
  dim(cmx_6) <- c(N_age, N_age)
  dim(cmx_7) <- c(N_age, N_age)
  dim(cmx_8) <- c(N_age, N_age)
  dim(cmx_9) <- c(N_age, N_age)
  
  dim(contact_matrix) <- c(N_age, N_age)
  
  dim(lambda) <- c(N_age,4)
  dim(lambda_1) <- N_age
  dim(lambda_2) <- N_age
  dim(lambda_3) <- N_age
  dim(lambda_4) <- N_age
  
  
  
  
})




