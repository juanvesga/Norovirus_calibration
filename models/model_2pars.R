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
    n_muM[i] -
    n_Mvacc[i]  
  
  
  update(S[]) <-
    S[i] + 
    n_bS[i] - 
    n_ageoS[i] + 
    n_ageiS[i] + 
    n_MS[i] + 
    sum(n_R_wane[i,])+
    sum(n_R1_wane[i,])+
    sum(n_R2_wane[i,])+
    sum(n_R3_wane[i,])+
    sum(n_R4_wane[i,])+
    sum(n_R12_wane[i,])+
    sum(n_R13_wane[i,])+
    sum(n_R14_wane[i,])+
    sum(n_R23_wane[i,])+
    sum(n_R24_wane[i,])+
    sum(n_R34_wane[i,])+
    sum(n_R4rd_wane[i,])-
    sum(n_S_E[i,]) -
    n_muS[i] -
    n_Svacc[i]
  
  
  dim(n_S_E)<-c(N_age,N_strain)
  # Infection # 1: age x 4 dimensions
  update(E[,]) <-
    E[i,j]- 
    n_ageoE[i,j] + 
    n_ageiE[i,j] + 
    n_S_E[i,j]  - 
    n_E_I[i,j] - 
    n_muE[i,j] -
    n_E_vacc[i,j]
  
  dim(n_E_I)<-c(N_age,N_strain)
  dim(n_E_vacc)<-c(N_age,N_strain)
  
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
    n_A_R[i,j] - 
    sum(n_R_tot[i,j,])-
    n_R_wane[i,j]-
    n_muR[i,j]-
    n_Rvacc[i,j]
  
  dim(n_R_wane)<-c(N_age,N_strain)
  
  ############## 2nd infection coming from GI3 (1)
  update(E1[,]) <-
    E1[i,j]- 
    n_ageoE1[i,j] + 
    n_ageiE1[i,j] +
    n_R_E1[i,j] - 
    n_E1_I1[i,j] - 
    n_muE1[i,j] - 
    n_E1_vacc[i,j] 
  
  
  update(I1[,]) <-   
    I1[i,j] - 
    n_ageoI1[i,j] + 
    n_ageiI1[i,j] + 
    n_E1_I1[i,j] -
    n_I1_A1[i,j] - 
    n_muI1[i,j] 
  
  
  update(A1[,]) <-  
    A1[i,j]-         
    n_ageoA1[i,j] + # age out
    n_ageiA1[i,j] + # age in
    n_I1_A1[i,j] +  # I to A 
    n_reinf_A1[i,j]+
    n_vacA1[i,j] +
    n_R1_A1[i,j] -  # Asymp re infection
    n_A1_R1[i,j] -  # Recovery A to R
    n_muA1[i,j]   
  
  
  update(R1[,]) <-  
    R1[i,j]- 
    n_ageoR1[i,j] + # age out
    n_ageiR1[i,j] + # age in
    n_A1_R1[i,j] -  # A to R recovery 
    sum(n_R1_tot[i,j,]) - # Third fresh infection
    n_R1_wane[i,j]-
    n_muR1[i,j] -  
    n_R1vacc[i,j] +
    n_vacR1[i,j]
  
  
  dim(n_R_E1)<-c(N_age,N_strain_1)
  dim(n_E1_I1)<-c(N_age,N_strain_1)
  dim(n_E1_vacc)<-c(N_age,N_strain_1)
  dim(n_I1_A1)<-c(N_age,N_strain_1)
  dim(n_A1_R1)<-c(N_age,N_strain_1)
  dim(n_R1_A1)<-c(N_age,N_strain_1)
  dim(n_R1_wane)<-c(N_age,N_strain_1)
  dim(n_reinf_A1)<-c(N_age,N_strain_1)
  
  ############### 2nd infection coming from other GI (2)
  update(E2[,]) <-
    E2[i,j]- 
    n_ageoE2[i,j] + 
    n_ageiE2[i,j] +
    n_R_E2[i,j] - 
    n_E2_I2[i,j] - 
    n_muE2[i,j] -
    n_E2_vacc[i,j]
  
  update(I2[,]) <-   
    I2[i,j] - 
    n_ageoI2[i,j] + 
    n_ageiI2[i,j] + 
    n_E2_I2[i,j] -
    n_I2_A2[i,j] - 
    n_muI2[i,j]  
  
  update(A2[,]) <-  
    A2[i,j]-         
    n_ageoA2[i,j] + # age out
    n_ageiA2[i,j] + # age in
    n_I2_A2[i,j] +  # I to A 
    n_reinf_A2[i,j]+
    n_R2_A2[i,j] -  # Asymp re infection
    n_A2_R2[i,j] -  # Recovery A to R
    n_muA2[i,j]    
  
  
  update(R2[,]) <-  
    R2[i,j]- 
    n_ageoR2[i,j] + # age out
    n_ageiR2[i,j] + # age in
    n_A2_R2[i,j] -  # A to R recovery 
    sum(n_R2_tot[i,j,]) - # Third fresh infection
    n_R2_wane[i,j]-
    n_muR2[i,j]    -  
    n_R2vacc[i,j] 
  
  dim(n_R_E2)<-c(N_age,N_strain_1)
  dim(n_E2_I2)<-c(N_age,N_strain_1)
  dim(n_E2_vacc)<-c(N_age,N_strain_1)
  dim(n_I2_A2)<-c(N_age,N_strain_1)
  dim(n_A2_R2)<-c(N_age,N_strain_1)
  dim(n_R2_A2)<-c(N_age,N_strain_1)
  dim(n_R2_wane)<-c(N_age,N_strain_1)
  dim(n_reinf_A2)<-c(N_age,N_strain_1)
  
  
  ############### 2nd infection coming from GIi4 (3)
  update(E3[,]) <-
    E3[i,j]- 
    n_ageoE3[i,j] + 
    n_ageiE3[i,j] +
    n_R_E3[i,j] - 
    n_E3_I3[i,j] - 
    n_muE3[i,j] -
    n_E3_vacc[i,j]
  
  update(I3[,]) <-   
    I3[i,j] - 
    n_ageoI3[i,j] + 
    n_ageiI3[i,j] + 
    n_E3_I3[i,j] -
    n_I3_A3[i,j] - 
    n_muI3[i,j] 
  
  update(A3[,]) <-  
    A3[i,j]-         
    n_ageoA3[i,j] + # age out
    n_ageiA3[i,j] + # age in
    n_I3_A3[i,j] +  # I to A 
    n_reinf_A3[i,j]+
    n_vacA3[i,j] +
    n_R3_A3[i,j] -  # Asymp re infection
    n_A3_R3[i,j] -  # Recovery A to R
    n_muA3[i,j]     # BAck. mortality
  
  
  update(R3[,]) <-  
    R3[i,j]- 
    n_ageoR3[i,j] + # age out
    n_ageiR3[i,j] + # age in
    n_A3_R3[i,j] -  # A to R recovery 
    sum(n_R3_tot[i,j,]) - # Third fresh infection
    n_R3_wane[i,j]-
    n_muR3[i,j] -  
    n_R3vacc[i,j]+
    n_vacR3[i,j]
  
  dim(n_R_E3)<-c(N_age,N_strain_1)
  dim(n_E3_I3)<-c(N_age,N_strain_1)
  dim(n_E3_vacc)<-c(N_age,N_strain_1)
  dim(n_I3_A3)<-c(N_age,N_strain_1)
  dim(n_A3_R3)<-c(N_age,N_strain_1)
  dim(n_R3_A3)<-c(N_age,N_strain_1)
  dim(n_R3_wane)<-c(N_age,N_strain_1)
  dim(n_reinf_A3)<-c(N_age,N_strain_1)
  
  ############### 2nd infection coming from other GII (4)
  update(E4[,]) <-
    E4[i,j]- 
    n_ageoE4[i,j] + 
    n_ageiE4[i,j] +
    n_R_E4[i,j] - 
    n_E4_I4[i,j] - 
    n_muE4[i,j] -
    n_E4_vacc[i,j]
  
  update(I4[,]) <-   
    I4[i,j] - 
    n_ageoI4[i,j] + 
    n_ageiI4[i,j] + 
    n_E4_I4[i,j] -
    n_I4_A4[i,j] - 
    n_muI4[i,j] 
  
  update(A4[,]) <-  
    A4[i,j]-         
    n_ageoA4[i,j] + # age out
    n_ageiA4[i,j] + # age in
    n_I4_A4[i,j] +  # I to A 
    n_reinf_A4[i,j]+
    n_R4_A4[i,j] -  # Asymp re infection
    n_A4_R4[i,j] -  # Recovery A to R
    n_muA4[i,j]     # BAck. mortality
  
  
  update(R4[,]) <-  
    R4[i,j]- 
    n_ageoR4[i,j] + # age out
    n_ageiR4[i,j] + # age in
    n_A4_R4[i,j] -  # A to R recovery 
    sum(n_R4_tot[i,j,]) - # Third fresh infection
    n_R4_wane[i,j]-
    n_muR4[i,j] -  
    n_R4vacc[i,j]  
  
  dim(n_R_E4)<-c(N_age,N_strain_1)
  dim(n_E4_I4)<-c(N_age,N_strain_1)
  dim(n_E4_vacc)<-c(N_age,N_strain_1)
  dim(n_I4_A4)<-c(N_age,N_strain_1)
  dim(n_A4_R4)<-c(N_age,N_strain_1)
  dim(n_R4_A4)<-c(N_age,N_strain_1)
  dim(n_R4_wane)<-c(N_age,N_strain_1)
  dim(n_reinf_A4)<-c(N_age,N_strain_1)
  
  #########HERE
  
  
  ############## 3rd infection coming from GI3 (1) & GI (2)
  update(E12[,]) <-
    E12[i,j]- 
    n_ageoE12[i,j] + # age out
    n_ageiE12[i,j] + # age in
    n_R1_E12[i,j] +
    n_R2_E12[i,j] - # Symp infection from R second infection (from 1 and 2 to 3 & 4)
    n_E12_I12[i,j] - 
    n_muE12[i,j] -
    n_E12_vacc[i,j]
  
  update(I12[,]) <-   
    I12[i,j] - 
    n_ageoI12[i,j] + 
    n_ageiI12[i,j] + 
    n_E12_I12[i,j] -
    n_I12_A12[i,j] - 
    n_muI12[i,j] 
  
  update(A12[,]) <-  
    A12[i,j]-         
    n_ageoA12[i,j] + # age out
    n_ageiA12[i,j] + # age in
    n_I12_A12[i,j] +  # I to A 
    n_reinf_A12[i,j]+
    n_vacA12[i,j] +
    n_R12_A12[i,j] -  # Asymp re infection
    n_A12_R12[i,j] -  # Recovery A to R
    n_muA12[i,j]     # BAck. mortality
  
  
  update(R12[,]) <-  
    R12[i,j]- 
    n_ageoR12[i,j] + # age out
    n_ageiR12[i,j] + # age in
    n_A12_R12[i,j] -  # A to R recovery 
    sum(n_R12_tot[i,j,]) - # fourth fresh infection
    n_R12_wane[i,j]-
    n_muR12[i,j] -  
    n_R12vacc[i,j] +
    n_vacR12[i,j]
  
  dim(n_R1_E12)<-c(N_age,N_strain_2)
  dim(n_R2_E12)<-c(N_age,N_strain_2)
  dim(n_E12_I12)<-c(N_age,N_strain_2)
  dim(n_E12_vacc)<-c(N_age,N_strain_2)
  dim(n_I12_A12)<-c(N_age,N_strain_2)
  dim(n_A12_R12)<-c(N_age,N_strain_2)
  dim(n_R12_A12)<-c(N_age,N_strain_2)
  dim(n_R12_wane)<-c(N_age,N_strain_2)
  dim(n_reinf_A12)<-c(N_age,N_strain_2)
  
  
  ############## 3rd infection coming from GI3 (1) & GII4 (3)
  update(E13[,]) <-
    E13[i,j]- 
    n_ageoE13[i,j] + # age out
    n_ageiE13[i,j] + # age in
    n_R1_E13[i,j]+
    n_R3_E13[i,j]- # Symp infection from R second infection (from 1 and 3 to 2 & 4)
    n_E13_I13[i,j] - 
    n_muE13[i,j] -
    n_E13_vacc[i,j]
  
  update(I13[,]) <-   
    I13[i,j] - 
    n_ageoI13[i,j] + 
    n_ageiI13[i,j] + 
    n_E13_I13[i,j] -
    n_I13_A13[i,j] - 
    n_muI13[i,j] 
  
  update(A13[,]) <-  
    A13[i,j]-         
    n_ageoA13[i,j] + # age out
    n_ageiA13[i,j] + # age in
    n_I13_A13[i,j] +  # I to A 
    n_reinf_A13[i,j]+
    n_vacA13[i,j] +
    n_R13_A13[i,j] -  # Asymp re infection
    n_A13_R13[i,j] -  # Recovery A to R
    n_muA13[i,j]     # BAck. mortality
  
  
  update(R13[,]) <-  
    R13[i,j]- 
    n_ageoR13[i,j] + # age out
    n_ageiR13[i,j] + # age in
    n_A13_R13[i,j] -  # A to R recovery 
    sum(n_R13_tot[i,j,]) - # fourth fresh infection
    n_R13_wane[i,j]-
    n_muR13[i,j] -  
    n_R13vacc[i,j]  +
    n_vacR13[i,j]
  
  
  dim(n_R1_E13)<-c(N_age,N_strain_2)
  dim(n_R3_E13)<-c(N_age,N_strain_2)
  dim(n_E13_I13)<-c(N_age,N_strain_2)
  dim(n_E13_vacc)<-c(N_age,N_strain_2)
  dim(n_I13_A13)<-c(N_age,N_strain_2)
  dim(n_A13_R13)<-c(N_age,N_strain_2)
  dim(n_R13_A13)<-c(N_age,N_strain_2)
  dim(n_R13_wane)<-c(N_age,N_strain_2)
  dim(n_reinf_A13)<-c(N_age,N_strain_2)
  
  ############## 3rd infection coming from GI3 (1) & GII (4)
  update(E14[,]) <-
    E14[i,j]- 
    n_ageoE14[i,j] + # age out
    n_ageiE14[i,j] + # age in
    n_R1_E14[i,j]+
    n_R4_E14[i,j]- # Symp infection from R second infection (from 1 and 4 to 2 & 3)
    n_E14_I14[i,j] - 
    n_muE14[i,j] -
    n_E14_vacc[i,j]
  
  update(I14[,]) <-   
    I14[i,j] - 
    n_ageoI14[i,j] + 
    n_ageiI14[i,j] + 
    n_E14_I14[i,j] -
    n_I14_A14[i,j] - 
    n_muI14[i,j] 
  
  update(A14[,]) <-  
    A14[i,j]-         
    n_ageoA14[i,j] + # age out
    n_ageiA14[i,j] + # age in
    n_I14_A14[i,j] +  # I to A 
    n_reinf_A14[i,j]+
    n_vacA14[i,j] +
    n_R14_A14[i,j] -  # Asymp re infection
    n_A14_R14[i,j] -  # Recovery A to R
    n_muA14[i,j]     # BAck. mortality
  
  
  update(R14[,]) <-  
    R14[i,j]- 
    n_ageoR14[i,j] + # age out
    n_ageiR14[i,j] + # age in
    n_A14_R14[i,j] -  # A to R recovery 
    sum(n_R14_tot[i,j,]) - # fourth fresh infection
    n_R14_wane[i,j]-
    n_muR14[i,j] -  
    n_R14vacc[i,j] +
    n_vacR14[i,j]
  
  
  dim(n_R1_E14)<-c(N_age,N_strain_2)
  dim(n_R4_E14)<-c(N_age,N_strain_2)
  dim(n_E14_I14)<-c(N_age,N_strain_2)
  dim(n_E14_vacc)<-c(N_age,N_strain_2)
  dim(n_I14_A14)<-c(N_age,N_strain_2)
  dim(n_A14_R14)<-c(N_age,N_strain_2)
  dim(n_R14_A14)<-c(N_age,N_strain_2)
  dim(n_R14_wane)<-c(N_age,N_strain_2)
  dim(n_reinf_A14)<-c(N_age,N_strain_2)
  
  ############## 3rd infection coming from GI (2) & GII4 (3)
  update(E23[,]) <-
    E23[i,j]- 
    n_ageoE23[i,j] + # age out
    n_ageiE23[i,j] + # age in
    n_R2_E23[i,j] +
    n_R3_E23[i,j]- # Symp infection from R second infection (from 2 and 3 to 1 & 4)
    n_E23_I23[i,j] - 
    n_muE23[i,j] -
    n_E23_vacc[i,j]
  
  update(I23[,]) <-   
    I23[i,j] - 
    n_ageoI23[i,j] + 
    n_ageiI23[i,j] + 
    n_E23_I23[i,j] -
    n_I23_A23[i,j] - 
    n_muI23[i,j] 
  
  update(A23[,]) <-  
    A23[i,j]-         
    n_ageoA23[i,j] + # age out
    n_ageiA23[i,j] + # age in
    n_I23_A23[i,j] +  # I to A 
    n_reinf_A23[i,j]+
    n_vacA23[i,j] +
    n_R23_A23[i,j] -  # Asymp re infection
    n_A23_R23[i,j] -  # Recovery A to R
    n_muA23[i,j]     # BAck. mortality
  
  update(R23[,]) <-  
    R23[i,j]- 
    n_ageoR23[i,j] + # age out
    n_ageiR23[i,j] + # age in
    n_A23_R23[i,j] -  # A to R recovery 
    sum(n_R23_tot[i,j,]) - # fourth fresh infection
    n_R23_wane[i,j]-
    n_muR23[i,j] -  
    n_R23vacc[i,j] +
    n_vacR23[i,j]
  
  dim(n_R2_E23)<-c(N_age,N_strain_2)
  dim(n_R3_E23)<-c(N_age,N_strain_2)
  dim(n_E23_I23)<-c(N_age,N_strain_2)
  dim(n_E23_vacc)<-c(N_age,N_strain_2)
  dim(n_I23_A23)<-c(N_age,N_strain_2)
  dim(n_A23_R23)<-c(N_age,N_strain_2)
  dim(n_R23_A23)<-c(N_age,N_strain_2)
  dim(n_R23_wane)<-c(N_age,N_strain_2)
  dim(n_reinf_A23)<-c(N_age,N_strain_2)
  
  ############## 3rd infection coming from GI (2) & GII (4)
  update(E24[,]) <-
    E24[i,j]- 
    n_ageoE24[i,j] + # age out
    n_ageiE24[i,j] + # age in
    n_R2_E24[i,j]+
    n_R4_E24[i,j]- # Symp infection from R second infection (from 2 and 4 to 1 & 3)
    n_E24_I24[i,j] - 
    n_muE24[i,j] -
    n_E24_vacc[i,j]
  
  update(I24[,]) <-   
    I24[i,j] - 
    n_ageoI24[i,j] + 
    n_ageiI24[i,j] + 
    n_E24_I24[i,j] -
    n_I24_A24[i,j] - 
    n_muI24[i,j] 
  
  update(A24[,]) <-  
    A24[i,j]-         
    n_ageoA24[i,j] + # age out
    n_ageiA24[i,j] + # age in
    n_I24_A24[i,j] +  # I to A 
    n_reinf_A24[i,j]+
    n_R24_A24[i,j] -  # Asymp re infection
    n_A24_R24[i,j] -  # Recovery A to R
    n_muA24[i,j]     # BAck. mortality
  
  
  update(R24[,]) <-  
    R24[i,j]- 
    n_ageoR24[i,j] + # age out
    n_ageiR24[i,j] + # age in
    n_A24_R24[i,j] -  # A to R recovery 
    sum(n_R24_tot[i,j,]) - # fourth fresh infection
    n_R24_wane[i,j]-
    n_muR24[i,j] -  
    n_R24vacc[i,j]
  
  dim(n_R2_E24)<-c(N_age,N_strain_2)
  dim(n_R4_E24)<-c(N_age,N_strain_2)
  dim(n_E24_I24)<-c(N_age,N_strain_2)
  dim(n_E24_vacc)<-c(N_age,N_strain_2)
  dim(n_I24_A24)<-c(N_age,N_strain_2)
  dim(n_A24_R24)<-c(N_age,N_strain_2)
  dim(n_R24_A24)<-c(N_age,N_strain_2)
  dim(n_R24_wane)<-c(N_age,N_strain_2)
  dim(n_reinf_A24)<-c(N_age,N_strain_2)
  
  ############## 3rd infection coming from GII4 (3) & GII (4)
  update(E34[,]) <-
    E34[i,j]- 
    n_ageoE34[i,j] + # age out
    n_ageiE34[i,j] + # age in
    n_R3_E34[i,j]+
    n_R4_E34[i,j]- # Symp infection from R second infection (from 3 and 4 to 1 & 2)
    n_E34_I34[i,j] - 
    n_muE34[i,j] -
    n_E34_vacc[i,j] 
  
  update(I34[,]) <-   
    I34[i,j] - 
    n_ageoI34[i,j] + 
    n_ageiI34[i,j] + 
    n_E34_I34[i,j] -
    n_I34_A34[i,j] - 
    n_muI34[i,j] 
  
  update(A34[,]) <-  
    A34[i,j]-         
    n_ageoA34[i,j] + # age out
    n_ageiA34[i,j] + # age in
    n_I34_A34[i,j] +  # I to A 
    n_reinf_A34[i,j]+
    n_vacA34[i,j] +
    n_R34_A34[i,j] -  # Asymp re infection
    n_A34_R34[i,j] -  # Recovery A to R
    n_muA34[i,j]     # Back. mortality
  
  
  update(R34[,]) <-  
    R34[i,j]- 
    n_ageoR34[i,j] + # age out
    n_ageiR34[i,j] + # age in
    n_A34_R34[i,j] -  # A to R recovery 
    sum(n_R34_tot[i,j,]) - # fourth fresh infection
    n_R34_wane[i,j]-
    n_muR34[i,j]-  
    n_R34vacc[i,j] +
    n_vacR34[i,j]
  
  dim(n_R3_E34)<-c(N_age,N_strain_2)
  dim(n_R4_E34)<-c(N_age,N_strain_2)
  dim(n_E34_I34)<-c(N_age,N_strain_2)
  dim(n_E34_vacc)<-c(N_age,N_strain_2)
  dim(n_I34_A34)<-c(N_age,N_strain_2)
  dim(n_A34_R34)<-c(N_age,N_strain_2)
  dim(n_R34_A34)<-c(N_age,N_strain_2)
  dim(n_R34_wane)<-c(N_age,N_strain_2)
  dim(n_reinf_A34)<-c(N_age,N_strain_2)
  
  
  
  ############## 4rd infection 
  update(E4rd[,]) <-
    E4rd[i,j]- 
    n_ageoE4rd[i,j] + # age out
    n_ageiE4rd[i,j] + # age in
    n_R12_E4rd[i,j] +
    n_R13_E4rd[i,j] +
    n_R14_E4rd[i,j] +
    n_R23_E4rd[i,j] +
    n_R24_E4rd[i,j] +
    n_R34_E4rd[i,j] - # Symp infection from R second infection (from 3 and 4 to 1 & 2)
    n_E4rd_I4rd[i,j] - 
    n_muE4rd[i,j] -
    n_E4rd_vacc[i,j]
  
  update(I4rd[,]) <-   
    I4rd[i,j] - 
    n_ageoI4rd[i,j] + 
    n_ageiI4rd[i,j] + 
    n_E4rd_I4rd[i,j] -
    n_I4rd_A4rd[i,j] - 
    n_muI4rd[i,j] 
  
  update(A4rd[,]) <-  
    A4rd[i,j]-         
    n_ageoA4rd[i,j] + # age out
    n_ageiA4rd[i,j] + # age in
    n_I4rd_A4rd[i,j] +  # I to A 
     n_vacA1234[i,j] +
    n_R4rd_A4rd[i,j] -  # Asymp re infection
    n_A4rd_R4rd[i,j] -  # Recovery A to R
    n_muA4rd[i,j]     # Back. mortality
  
  
  update(R4rd[,]) <-  
    R4rd[i,j]- 
    n_ageoR4rd[i,j] + # age out
    n_ageiR4rd[i,j] + # age in
    n_A4rd_R4rd[i,j] -  # A to R recovery 
    sum(n_R4rd_tot[i,j,]) -  # All reinf
    n_R4rd_wane[i,j]-
    n_muR4rd[i,j]-  
    n_R4rdvacc[i,j]+
    n_vacR1234[i,j]
  
  
  dim(n_R12_E4rd)<-c(N_age,N_strain)
  dim(n_R13_E4rd)<-c(N_age,N_strain)
  dim(n_R14_E4rd)<-c(N_age,N_strain)
  dim(n_R23_E4rd)<-c(N_age,N_strain)
  dim(n_R24_E4rd)<-c(N_age,N_strain)
  dim(n_R34_E4rd)<-c(N_age,N_strain)
  
  dim(n_E4rd_I4rd)<-c(N_age,N_strain)
  dim(n_E4rd_vacc)<-c(N_age,N_strain)
  dim(n_I4rd_A4rd)<-c(N_age,N_strain)
  dim(n_A4rd_R4rd)<-c(N_age,N_strain)
  dim(n_R4rd_A4rd)<-c(N_age,N_strain)
  dim(n_R4rd_wane)<-c(N_age,N_strain)
  
  
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
         sum(n_E_I[1:4,1]) * reported_fac[i] +
         sum(n_E2_I2[1:4,1]) * reported_fac[i] +
         sum(n_E3_I3[1:4,1]) * reported_fac[i] +
         sum(n_E4_I4[1:4,1]) * reported_fac[i] +
         sum(n_E23_I23[1:4,1]) * reported_fac[i] +
         sum(n_E24_I24[1:4,1]) * reported_fac[i] +
         sum(n_E34_I34[1:4,1]) * reported_fac[i] +
         sum(n_E4rd_I4rd[1:4,1]) * reported_fac[i]  
      ) else
        if(i==2) 
          (infections_day_gi3[i]+
             sum(n_E_I[5:8,1]) * reported_fac[i] +
             sum(n_E2_I2[5:8,1]) * reported_fac[i] +
             sum(n_E3_I3[5:8,1]) * reported_fac[i] +
             sum(n_E4_I4[5:8,1]) * reported_fac[i] +
             sum(n_E23_I23[5:8,1]) * reported_fac[i] +
             sum(n_E24_I24[5:8,1]) * reported_fac[i] +
             sum(n_E34_I34[5:8,1]) * reported_fac[i] +
             sum(n_E4rd_I4rd[5:8,1]) * reported_fac[i]  
          ) else
            if(i==3) 
              (infections_day_gi3[i]+
                 sum(n_E_I[9:13,1]) * reported_fac[i] +
                 sum(n_E2_I2[9:13,1]) * reported_fac[i] +
                 sum(n_E3_I3[9:13,1]) * reported_fac[i] +
                 sum(n_E4_I4[9:13,1]) * reported_fac[i] +
                 sum(n_E23_I23[9:13,1]) * reported_fac[i] +
                 sum(n_E24_I24[9:13,1]) * reported_fac[i] +
                 sum(n_E34_I34[9:13,1]) * reported_fac[i] +
                 sum(n_E4rd_I4rd[9:13,1]) * reported_fac[i]  
              ) else
                (infections_day_gi3[i]+
                   sum(n_E_I[14,1]) * reported_fac[i] +
                   sum(n_E2_I2[14,1]) * reported_fac[i] +
                   sum(n_E3_I3[14,1]) * reported_fac[i] +
                   sum(n_E4_I4[14,1]) * reported_fac[i] +
                   sum(n_E23_I23[14,1]) * reported_fac[i] +
                   sum(n_E24_I24[14,1]) * reported_fac[i] +
                   sum(n_E34_I34[14,1]) * reported_fac[i] +
                   sum(n_E4rd_I4rd[14,1]) * reported_fac[i])
  
  
  update(infections_day_gi[]) <- 
    if (i==1)
      (infections_day_gi[i]+
         sum(n_E_I[1:4,2]) * reported_fac[i] +
         sum(n_E1_I1[1:4,1]) * reported_fac[i] +
         sum(n_E3_I3[1:4,2]) * reported_fac[i] +
         sum(n_E4_I4[1:4,2]) * reported_fac[i] +
         sum(n_E13_I13[1:4,1]) * reported_fac[i] +
         sum(n_E14_I14[1:4,1]) * reported_fac[i] +
         sum(n_E34_I34[1:4,2]) * reported_fac[i] +
         sum(n_E4rd_I4rd[1:4,2]) * reported_fac[i]  
      ) else
        if(i==2) 
          (infections_day_gi[i]+
             sum(n_E_I[5:8,2]) * reported_fac[i] +
             sum(n_E1_I1[5:8,1]) * reported_fac[i] +
             sum(n_E3_I3[5:8,2]) * reported_fac[i] +
             sum(n_E4_I4[5:8,2]) * reported_fac[i] +
             sum(n_E13_I13[5:8,1]) * reported_fac[i] +
             sum(n_E14_I14[5:8,1]) * reported_fac[i] +
             sum(n_E34_I34[5:8,2]) * reported_fac[i] +
             sum(n_E4rd_I4rd[5:8,2]) * reported_fac[i]  
          ) else
            if(i==3) 
              (infections_day_gi[i]+
                 sum(n_E_I[9:13,2]) * reported_fac[i] +
                 sum(n_E1_I1[9:13,1]) * reported_fac[i] +
                 sum(n_E3_I3[9:13,2]) * reported_fac[i] +
                 sum(n_E4_I4[9:13,2]) * reported_fac[i] +
                 sum(n_E13_I13[9:13,1]) * reported_fac[i] +
                 sum(n_E14_I14[9:13,1]) * reported_fac[i] +
                 sum(n_E34_I34[9:13,2]) * reported_fac[i] +
                 sum(n_E4rd_I4rd[9:13,2]) * reported_fac[i] 
              ) else
                (infections_day_gi[i]+
                   sum(n_E_I[14,2]) * reported_fac[i] +
                   sum(n_E1_I1[14,1]) * reported_fac[i] +
                   sum(n_E3_I3[14,2]) * reported_fac[i] +
                   sum(n_E4_I4[14,2]) * reported_fac[i] +
                   sum(n_E13_I13[14,1]) * reported_fac[i] +
                   sum(n_E14_I14[14,1]) * reported_fac[i] +
                   sum(n_E34_I34[14,2]) * reported_fac[i] +
                   sum(n_E4rd_I4rd[14,2]) * reported_fac[i] )
  
  
  update(infections_day_gii4[]) <- 
    if (i==1)
      (infections_day_gii4[i]+
         sum(n_E_I[1:4,3]) * reported_fac[i] +
         sum(n_E1_I1[1:4,2]) * reported_fac[i] +
         sum(n_E2_I2[1:4,2]) * reported_fac[i] +
         sum(n_E4_I4[1:4,3]) * reported_fac[i] +
         sum(n_E12_I12[1:4,1]) * reported_fac[i] +
         sum(n_E14_I14[1:4,2]) * reported_fac[i] +
         sum(n_E24_I24[1:4,2]) * reported_fac[i] +
         sum(n_E4rd_I4rd[1:4,3]) * reported_fac[i]  
      ) else
        if(i==2) 
          (infections_day_gii4[i]+
             sum(n_E_I[5:8,3]) * reported_fac[i] +
             sum(n_E1_I1[5:8,2]) * reported_fac[i] +
             sum(n_E2_I2[5:8,2]) * reported_fac[i] +
             sum(n_E4_I4[5:8,3]) * reported_fac[i] +
             sum(n_E12_I12[5:8,1]) * reported_fac[i] +
             sum(n_E14_I14[5:8,2]) * reported_fac[i] +
             sum(n_E24_I24[5:8,2]) * reported_fac[i] +
             sum(n_E4rd_I4rd[5:8,3]) * reported_fac[i]  
          ) else
            if(i==3) 
              (infections_day_gii4[i]+
                 sum(n_E_I[9:13,3]) * reported_fac[i] +
                 sum(n_E1_I1[9:13,2]) * reported_fac[i] +
                 sum(n_E2_I2[9:13,2]) * reported_fac[i] +
                 sum(n_E4_I4[9:13,3]) * reported_fac[i] +
                 sum(n_E12_I12[9:13,1]) * reported_fac[i] +
                 sum(n_E14_I14[9:13,2]) * reported_fac[i] +
                 sum(n_E24_I24[9:13,2]) * reported_fac[i] +
                 sum(n_E4rd_I4rd[9:13,3]) * reported_fac[i]  
              ) else
                (infections_day_gii4[i]+
                   sum(n_E_I[14,3]) * reported_fac[i] +
                   sum(n_E1_I1[14,2]) * reported_fac[i] +
                   sum(n_E2_I2[14,2]) * reported_fac[i] +
                   sum(n_E4_I4[14,3]) * reported_fac[i] +
                   sum(n_E12_I12[14,1]) * reported_fac[i] +
                   sum(n_E14_I14[14,2]) * reported_fac[i] +
                   sum(n_E24_I24[14,2]) * reported_fac[i] +
                   sum(n_E4rd_I4rd[14,3]) * reported_fac[i])
  
  
  
  update(infections_day_gii[]) <- 
    if (i==1)
      (infections_day_gii[i]+
         sum(n_E_I[1:4,4]) * reported_fac[i] +
         sum(n_E1_I1[1:4,3]) * reported_fac[i] +
         sum(n_E2_I2[1:4,3]) * reported_fac[i] +
         sum(n_E3_I3[1:4,3]) * reported_fac[i] +
         sum(n_E12_I12[1:4,2]) * reported_fac[i] +
         sum(n_E13_I13[1:4,2]) * reported_fac[i] +
         sum(n_E23_I23[1:4,2]) * reported_fac[i] +
         sum(n_E4rd_I4rd[1:4,4]) * reported_fac[i]  
      ) else
        if(i==2) 
          (infections_day_gii[i]+
             sum(n_E_I[5:8,4]) * reported_fac[i] +
             sum(n_E1_I1[5:8,3]) * reported_fac[i] +
             sum(n_E2_I2[5:8,3]) * reported_fac[i] +
             sum(n_E3_I3[5:8,3]) * reported_fac[i] +
             sum(n_E12_I12[5:8,2]) * reported_fac[i] +
             sum(n_E13_I13[5:8,2]) * reported_fac[i] +
             sum(n_E23_I23[5:8,2]) * reported_fac[i] +
             sum(n_E4rd_I4rd[5:8,4]) * reported_fac[i]    
          ) else
            if(i==3) 
              (infections_day_gii[i]+
                 sum(n_E_I[9:13,4]) * reported_fac[i] +
                 sum(n_E1_I1[9:13,3]) * reported_fac[i] +
                 sum(n_E2_I2[9:13,3]) * reported_fac[i] +
                 sum(n_E3_I3[9:13,3]) * reported_fac[i] +
                 sum(n_E12_I12[9:13,2]) * reported_fac[i] +
                 sum(n_E13_I13[9:13,2]) * reported_fac[i] +
                 sum(n_E23_I23[9:13,2]) * reported_fac[i] +
                 sum(n_E4rd_I4rd[9:13,4]) * reported_fac[i]   
              ) else
                (infections_day_gii[i]+
                   sum(n_E_I[14,4]) * reported_fac[i] +
                   sum(n_E1_I1[14,3]) * reported_fac[i] +
                   sum(n_E2_I2[14,3]) * reported_fac[i] +
                   sum(n_E3_I3[14,3]) * reported_fac[i] +
                   sum(n_E12_I12[14,2]) * reported_fac[i] +
                   sum(n_E13_I13[14,2]) * reported_fac[i] +
                   sum(n_E23_I23[14,2]) * reported_fac[i] +
                   sum(n_E4rd_I4rd[14,4]) * reported_fac[i]  )
  
  
  # Weekly reported cases (match sgss)
  
  update(reported_wk[]) <- reported_wk[i] + 
    (sum(n_E_I[i,]) * reported_fac_long[i] +
       sum(n_E1_I1[i,]) * reported_fac_long[i] +
       sum(n_E2_I2[i,]) * reported_fac_long[i] +
       sum(n_E3_I3[i,]) * reported_fac_long[i] +
       sum(n_E4_I4[i,]) * reported_fac_long[i] +
       sum(n_E12_I12[i,]) * reported_fac_long[i] +
       sum(n_E13_I13[i,]) * reported_fac_long[i] +
       sum(n_E14_I14[i,]) * reported_fac_long[i] +
       sum(n_E23_I23[i,]) * reported_fac_long[i] +
       sum(n_E24_I24[i,]) * reported_fac_long[i] +
       sum(n_E34_I34[i,]) * reported_fac_long[i] +
       sum(n_E4rd_I4rd[i,]) * reported_fac_long[i])  
  
  update(death_reported_wk_all[]) <- death_reported_wk_all[i] +
    (reported_wk[i] * p_cfr[i])
  
  update(death_reported_wk) <- death_reported_wk +
    sum(death_reported_wk_all)
  
  update(reported_wk_gi3) <- reported_wk_gi3 + 
    (sum(n_E_I[1:4,1])* reported_geno[1] +
       sum(n_E2_I2[1:4,1])* reported_geno[1] +
       sum(n_E3_I3[1:4,1])* reported_geno[1] +
       sum(n_E4_I4[1:4,1])* reported_geno[1] +
       sum(n_E23_I23[1:4,1])* reported_geno[1] +
       sum(n_E24_I24[1:4,1])* reported_geno[1] +
       sum(n_E34_I34[1:4,1])* reported_geno[1] +
       sum(n_E4rd_I4rd[1:4,1])* reported_geno[1] + 
       sum(n_E_I[5:8,1])* reported_geno[2] +
       sum(n_E2_I2[5:8,1])* reported_geno[2] +
       sum(n_E3_I3[5:8,1])* reported_geno[2] +
       sum(n_E4_I4[5:8,1])* reported_geno[2] +
       sum(n_E23_I23[5:8,1])* reported_geno[2] +
       sum(n_E24_I24[5:8,1])* reported_geno[2] +
       sum(n_E34_I34[5:8,1])* reported_geno[2] +
       sum(n_E4rd_I4rd[5:8,1])* reported_geno[2] +
       sum(n_E_I[9:13,1])* reported_geno[3] +
       sum(n_E2_I2[9:13,1])* reported_geno[3] +
       sum(n_E3_I3[9:13,1])* reported_geno[3] +
       sum(n_E4_I4[9:13,1])* reported_geno[3] +
       sum(n_E23_I23[9:13,1])* reported_geno[3] +
       sum(n_E24_I24[9:13,1])* reported_geno[3] +
       sum(n_E34_I34[9:13,1])* reported_geno[3] +
       sum(n_E4rd_I4rd[9:13,1])* reported_geno[3] +
       sum(n_E_I[14,1])* reported_geno[4] +
       sum(n_E2_I2[14,1])* reported_geno[4] +
       sum(n_E3_I3[14,1])* reported_geno[4] +
       sum(n_E4_I4[14,1])* reported_geno[4] +
       sum(n_E23_I23[14,1])* reported_geno[4] +
       sum(n_E24_I24[14,1])* reported_geno[4] +
       sum(n_E34_I34[14,1])* reported_geno[4] +
       sum(n_E4rd_I4rd[14,1])* reported_geno[4] ) 
  
  
  
  update(reported_wk_gi) <- reported_wk_gi + 
    (sum(n_E_I[1:4,2] )* reported_geno[1] +
       sum(n_E1_I1[1:4,1])* reported_geno[1]+
       sum(n_E3_I3[1:4,2])* reported_geno[1] +
       sum(n_E4_I4[1:4,2])* reported_geno[1] +
       sum(n_E13_I13[1:4,1])* reported_geno[1] +
       sum(n_E14_I14[1:4,1])* reported_geno[1] +
       sum(n_E34_I34[1:4,2])* reported_geno[1] +
       sum(n_E4rd_I4rd[1:4,2])* reported_geno[1] +
       sum(n_E_I[5:8,2] )* reported_geno[2] +
       sum(n_E1_I1[5:8,1])* reported_geno[2]+
       sum(n_E3_I3[5:8,2])* reported_geno[2] +
       sum(n_E4_I4[5:8,2])* reported_geno[2] +
       sum(n_E13_I13[5:8,1])* reported_geno[2] +
       sum(n_E14_I14[5:8,1])* reported_geno[2] +
       sum(n_E34_I34[5:8,2])* reported_geno[2] +
       sum(n_E4rd_I4rd[5:8,2])* reported_geno[2]+
       sum(n_E_I[9:13,2] )* reported_geno[3] +
       sum(n_E1_I1[9:13,1])* reported_geno[3]+
       sum(n_E3_I3[9:13,2])* reported_geno[3] +
       sum(n_E4_I4[9:13,2])* reported_geno[3] +
       sum(n_E13_I13[9:13,1])* reported_geno[3] +
       sum(n_E14_I14[9:13,1])* reported_geno[3] +
       sum(n_E34_I34[9:13,2])* reported_geno[3] +
       sum(n_E4rd_I4rd[9:13,2])* reported_geno[3]+
       sum(n_E_I[14,2] )* reported_geno[4] +
       sum(n_E1_I1[14,1])* reported_geno[4]+
       sum(n_E3_I3[14,2])* reported_geno[4] +
       sum(n_E4_I4[14,2])* reported_geno[4] +
       sum(n_E13_I13[14,1])* reported_geno[4] +
       sum(n_E14_I14[14,1])* reported_geno[4] +
       sum(n_E34_I34[14,2])* reported_geno[4] +
       sum(n_E4rd_I4rd[14,2])* reported_geno[4] ) 
  
  
  
  update(reported_wk_gii4) <- reported_wk_gii4 + 
    (sum(n_E_I[1:4,3])* reported_geno[1] +
       sum(n_E1_I1[1:4,2])* reported_geno[1] +
       sum(n_E2_I2[1:4,2])* reported_geno[1] +
       sum(n_E4_I4[1:4,3])* reported_geno[1] +
       sum(n_E12_I12[1:4,1])* reported_geno[1] +
       sum(n_E14_I14[1:4,2])* reported_geno[1] +
       sum(n_E24_I24[1:4,2])* reported_geno[1] +
       sum(n_E4rd_I4rd[1:4,3])* reported_geno[1] +
       sum(n_E_I[5:8,3])* reported_geno[2] +
       sum(n_E1_I1[5:8,2])* reported_geno[2] +
       sum(n_E2_I2[5:8,2])* reported_geno[2] +
       sum(n_E4_I4[5:8,3])* reported_geno[2] +
       sum(n_E12_I12[5:8,1])* reported_geno[2] +
       sum(n_E14_I14[5:8,2])* reported_geno[2] +
       sum(n_E24_I24[5:8,2])* reported_geno[2] +
       sum(n_E4rd_I4rd[5:8,3])* reported_geno[2]+
       sum(n_E_I[9:13,3])* reported_geno[3] +
       sum(n_E1_I1[9:13,2])* reported_geno[3] +
       sum(n_E2_I2[9:13,2])* reported_geno[3] +
       sum(n_E4_I4[9:13,3])* reported_geno[3] +
       sum(n_E12_I12[9:13,1])* reported_geno[3] +
       sum(n_E14_I14[9:13,2])* reported_geno[3] +
       sum(n_E24_I24[9:13,2])* reported_geno[3] +
       sum(n_E4rd_I4rd[9:13,3])* reported_geno[3]+
       sum(n_E_I[14,3])* reported_geno[4] +
       sum(n_E1_I1[14,2])* reported_geno[4] +
       sum(n_E2_I2[14,2])* reported_geno[4] +
       sum(n_E4_I4[14,3])* reported_geno[4] +
       sum(n_E12_I12[14,1])* reported_geno[4] +
       sum(n_E14_I14[14,2])* reported_geno[4] +
       sum(n_E24_I24[14,2])* reported_geno[4] +
       sum(n_E4rd_I4rd[14,3])* reported_geno[4]) 
  
  
  
  update(reported_wk_gii) <- reported_wk_gii + 
    (sum(n_E_I[1:4,4])* reported_geno[1] +
       sum(n_E1_I1[1:4,3])* reported_geno[1] +
       sum(n_E2_I2[1:4,3])* reported_geno[1] +
       sum(n_E3_I3[1:4,3])* reported_geno[1] +
       sum(n_E12_I12[1:4,2])* reported_geno[1] +
       sum(n_E13_I13[1:4,2])* reported_geno[1] +
       sum(n_E23_I23[1:4,2])* reported_geno[1] +
       sum(n_E4rd_I4rd[1:4,4])* reported_geno[1]  +
       sum(n_E_I[5:8,4])* reported_geno[2] +
       sum(n_E1_I1[5:8,3])* reported_geno[2] +
       sum(n_E2_I2[5:8,3])* reported_geno[2] +
       sum(n_E3_I3[5:8,3])* reported_geno[2] +
       sum(n_E12_I12[5:8,2])* reported_geno[2] +
       sum(n_E13_I13[5:8,2])* reported_geno[2] +
       sum(n_E23_I23[5:8,2])* reported_geno[2] +
       sum(n_E4rd_I4rd[5:8,4])* reported_geno[2]  +
       sum(n_E_I[9:13,4])* reported_geno[3] +
       sum(n_E1_I1[9:13,3])* reported_geno[3] +
       sum(n_E2_I2[9:13,3])* reported_geno[3] +
       sum(n_E3_I3[9:13,3])* reported_geno[3] +
       sum(n_E12_I12[9:13,2])* reported_geno[3] +
       sum(n_E13_I13[9:13,2])* reported_geno[3] +
       sum(n_E23_I23[9:13,2])* reported_geno[3] +
       sum(n_E4rd_I4rd[9:13,4])* reported_geno[3] +
       sum(n_E_I[14,4])* reported_geno[4] +
       sum(n_E1_I1[14,3])* reported_geno[4] +
       sum(n_E2_I2[14,3])* reported_geno[4] +
       sum(n_E3_I3[14,3])* reported_geno[4] +
       sum(n_E12_I12[14,2])* reported_geno[4] +
       sum(n_E13_I13[14,2])* reported_geno[4] +
       sum(n_E23_I23[14,2])* reported_geno[4] +
       sum(n_E4rd_I4rd[14,4])* reported_geno[4]) 
  
  
  
  
  
  
  update(inc_day_all_gi3[])<- 
    n_E_I[i,1]  +
    n_E2_I2[i,1]  +
    n_E3_I3[i,1] +
    n_E4_I4[i,1]  +
    n_E23_I23[i,1]  +
    n_E24_I24[i,1]  +
    n_E34_I34[i,1]  +
    n_E4rd_I4rd[i,1] 
  
  update(inc_day_all_gi[])<-
    n_E_I[i,2]  +
    n_E1_I1[i,1]  +
    n_E3_I3[i,2] +
    n_E4_I4[i,2]  +
    n_E13_I13[i,1]  +
    n_E14_I14[i,1] +
    n_E34_I34[i,2]  +
    n_E4rd_I4rd[i,2]   
  
  update(inc_day_all_gii4[])<- 
    n_E_I[i,3]  +
    n_E1_I1[i,2]  +
    n_E2_I2[i,2]  +
    n_E4_I4[i,3]  +
    n_E12_I12[i,1]  +
    n_E14_I14[i,2] +
    n_E24_I24[i,2]  +
    n_E4rd_I4rd[i,3] 
  
  update(inc_day_all_gii[])<- 
    n_E_I[i,4]  +
    n_E1_I1[i,3]  +
    n_E2_I2[i,3]  +
    n_E3_I3[i,3] +
    n_E12_I12[i,2]  +
    n_E13_I13[i,2]  +
    n_E23_I23[i,2]  +
    n_E4rd_I4rd[i,4] 
  
  
  # By 5 age groups
  
  update(inc_day_gi3[])<-
    (if(i==1) 
      sum(n_E_I[1,1])  +
       sum(n_E2_I2[1,1])  +
       sum(n_E3_I3[1,1]) +
       sum(n_E4_I4[1,1])  +
       sum(n_E23_I23[1,1])  +
       sum(n_E24_I24[1,1])  +
       sum(n_E34_I34[1,1])  +
       sum(n_E4rd_I4rd[1,1]) 
     else if(i==2)  
       sum(n_E_I[2:4,1])  +
       sum(n_E2_I2[2:4,1])  +
       sum(n_E3_I3[2:4,1]) +
       sum(n_E4_I4[2:4,1])  +
       sum(n_E23_I23[2:4,1])  +
       sum(n_E24_I24[2:4,1])  +
       sum(n_E34_I34[2:4,1])  +
       sum(n_E4rd_I4rd[2:4,1]) 
     else if(i==3)  
       sum(n_E_I[5:8,1])  +
       sum(n_E2_I2[5:8,1])  +
       sum(n_E3_I3[5:8,1]) +
       sum(n_E4_I4[5:8,1])  +
       sum(n_E23_I23[5:8,1])  +
       sum(n_E24_I24[5:8,1])  +
       sum(n_E34_I34[5:8,1])  +
       sum(n_E4rd_I4rd[5:8,1]) 
     else if (i==4) 
       sum(n_E_I[9:13,1])  +
       sum(n_E2_I2[9:13,1])  +
       sum(n_E3_I3[9:13,1]) +
       sum(n_E4_I4[9:13,1])  +
       sum(n_E23_I23[9:13,1])  +
       sum(n_E24_I24[9:13,1])  +
       sum(n_E34_I34[9:13,1])  +
       sum(n_E4rd_I4rd[9:13,1])   
     else
       sum(n_E_I[14,1])  +
       sum(n_E2_I2[14,1])  +
       sum(n_E3_I3[14,1]) +
       sum(n_E4_I4[14,1])  +
       sum(n_E23_I23[14,1])  +
       sum(n_E24_I24[14,1])  +
       sum(n_E34_I34[14,1])  +
       sum(n_E4rd_I4rd[14,1])  ) 
  
  
  update(inc_day_gi[])<-
    (if(i==1) 
      sum(n_E_I[1,2] ) +
       sum(n_E1_I1[1,1])  +
       sum(n_E3_I3[1,2]) +
       sum(n_E4_I4[1,2])  +
       sum(n_E13_I13[1,1])  +
       sum(n_E14_I14[1,1]) +
       sum(n_E34_I34[1,2])  +
       sum(n_E4rd_I4rd[1,2])  
     else if(i==2)  
       sum(n_E_I[2:4,2] ) +
       sum(n_E1_I1[2:4,1])  +
       sum(n_E3_I3[2:4,2]) +
       sum(n_E4_I4[2:4,2])  +
       sum(n_E13_I13[2:4,1])  +
       sum(n_E14_I14[2:4,1]) +
       sum(n_E34_I34[2:4,2])  +
       sum(n_E4rd_I4rd[2:4,2])
     else if(i==3)  
       sum(n_E_I[5:8,2] ) +
       sum(n_E1_I1[5:8,1])  +
       sum(n_E3_I3[5:8,2]) +
       sum(n_E4_I4[5:8,2])  +
       sum(n_E13_I13[5:8,1])  +
       sum(n_E14_I14[5:8,1]) +
       sum(n_E34_I34[5:8,2])  +
       sum(n_E4rd_I4rd[5:8,2])
     else if (i==4) 
       sum(n_E_I[9:13,2] ) +
       sum(n_E1_I1[9:13,1])  +
       sum(n_E3_I3[9:13,2]) +
       sum(n_E4_I4[9:13,2])  +
       sum(n_E13_I13[9:13,1])  +
       sum(n_E14_I14[9:13,1]) +
       sum(n_E34_I34[9:13,2])  +
       sum(n_E4rd_I4rd[9:13,2])
     else
       sum(n_E_I[14,2] ) +
       sum(n_E1_I1[14,1])  +
       sum(n_E3_I3[14,2]) +
       sum(n_E4_I4[14,2])  +
       sum(n_E13_I13[14,1])  +
       sum(n_E14_I14[14,1]) +
       sum(n_E34_I34[14,2])  +
       sum(n_E4rd_I4rd[14,2]) ) 
  
  update(inc_day_gii4[])<-
    (if(i==1) 
      sum(n_E_I[1,3])  +
       sum(n_E1_I1[1,2])  +
       sum(n_E2_I2[1,2])  +
       sum(n_E4_I4[1,3])  +
       sum(n_E12_I12[1,1])  +
       sum(n_E14_I14[1,2]) +
       sum(n_E24_I24[1,2])  +
       sum(n_E4rd_I4rd[1,3]) 
     else if(i==2)  
       sum(n_E_I[2:4,3])  +
       sum(n_E1_I1[2:4,2])  +
       sum(n_E2_I2[2:4,2])  +
       sum(n_E4_I4[2:4,3])  +
       sum(n_E12_I12[2:4,1])  +
       sum(n_E14_I14[2:4,2]) +
       sum(n_E24_I24[2:4,2])  +
       sum(n_E4rd_I4rd[2:4,3]) 
     else if(i==3)  
       sum(n_E_I[5:8,3])  +
       sum(n_E1_I1[5:8,2])  +
       sum(n_E2_I2[5:8,2])  +
       sum(n_E4_I4[5:8,3])  +
       sum(n_E12_I12[5:8,1])  +
       sum(n_E14_I14[5:8,2]) +
       sum(n_E24_I24[5:8,2])  +
       sum(n_E4rd_I4rd[5:8,3]) 
     else if (i==4) 
       sum(n_E_I[9:13,3])  +
       sum(n_E1_I1[9:13,2])  +
       sum(n_E2_I2[9:13,2])  +
       sum(n_E4_I4[9:13,3])  +
       sum(n_E12_I12[9:13,1])  +
       sum(n_E14_I14[9:13,2]) +
       sum(n_E24_I24[9:13,2])  +
       sum(n_E4rd_I4rd[9:13,3])    
     else
       sum(n_E_I[14,3])  +
       sum(n_E1_I1[14,2])  +
       sum(n_E2_I2[14,2])  +
       sum(n_E4_I4[14,3])  +
       sum(n_E12_I12[14,1])  +
       sum(n_E14_I14[14,2]) +
       sum(n_E24_I24[14,2])  +
       sum(n_E4rd_I4rd[14,3])  ) 
  
  
  update(inc_day_gii[])<-
    (if(i==1) 
      sum(n_E_I[1,4])  +
       sum(n_E1_I1[1,3])  +
       sum(n_E2_I2[1,3])  +
       sum(n_E3_I3[1,3]) +
       sum(n_E12_I12[1,2])  +
       sum(n_E13_I13[1,2])  +
       sum(n_E23_I23[1,2])  +
       sum(n_E4rd_I4rd[1,4]) 
     else if(i==2)  
       sum(n_E_I[2:4,4])  +
       sum(n_E1_I1[2:4,3])  +
       sum(n_E2_I2[2:4,3])  +
       sum(n_E3_I3[2:4,3]) +
       sum(n_E12_I12[2:4,2])  +
       sum(n_E13_I13[2:4,2])  +
       sum(n_E23_I23[2:4,2])  +
       sum(n_E4rd_I4rd[2:4,4])   
     else if(i==3)  
       sum(n_E_I[5:8,4])  +
       sum(n_E1_I1[5:8,3])  +
       sum(n_E2_I2[5:8,3])  +
       sum(n_E3_I3[5:8,3]) +
       sum(n_E12_I12[5:8,2])  +
       sum(n_E13_I13[5:8,2])  +
       sum(n_E23_I23[5:8,2])  +
       sum(n_E4rd_I4rd[5:8,4]) 
     else if (i==4) 
       sum(n_E_I[9:13,4])  +
       sum(n_E1_I1[9:13,3])  +
       sum(n_E2_I2[9:13,3])  +
       sum(n_E3_I3[9:13,3]) +
       sum(n_E12_I12[9:13,2])  +
       sum(n_E13_I13[9:13,2])  +
       sum(n_E23_I23[9:13,2])  +
       sum(n_E4rd_I4rd[9:13,4])   
     else
       sum(n_E_I[14,4])  +
       sum(n_E1_I1[14,3])  +
       sum(n_E2_I2[14,3])  +
       sum(n_E3_I3[14,3]) +
       sum(n_E12_I12[14,2])  +
       sum(n_E13_I13[14,2])  +
       sum(n_E23_I23[14,2])  +
       sum(n_E4rd_I4rd[14,4])  ) 
  
  
  
  dim(pop_forinc)<-5
  
  pop_forinc[1]<-N_byage[1]
  pop_forinc[2]<-sum(N_byage[2:4])
  pop_forinc[3]<-sum(N_byage[5:8])
  pop_forinc[4]<-sum(N_byage[9:13])
  pop_forinc[5]<-N_byage[14]
  
  dim(pop_all)<-N_age
  
  update(pop_all[])<-N_byage[i]
  
  
  
  
  update(inc_day_all[])<-1000*(inc_day_gi3[i]+inc_day_gi[i]+inc_day_gii4[i]+inc_day_gii[i])/pop_forinc[i]
  
  update(new_cases)<- sum(inc_day_gi3)+sum(inc_day_gi)+sum(inc_day_gii4)+sum(inc_day_gii)
  
  update(new_cases_week) <- new_cases_week + sum(inc_day_gi3)+sum(inc_day_gi)+sum(inc_day_gii4)+sum(inc_day_gii)
  
  update(new_hosp_adult_yr) <- new_hosp_adult_yr + 
    ( sum(inc_day_all_gii4[9:13]) +
    sum(inc_day_all_gii[9:13]) +
    sum(inc_day_all_gi3[9:13]) +
    sum(inc_day_all_gi[9:13]) ) * hosp_rate_adult
    
  
  update(new_hosp_elder_yr) <- new_hosp_elder_yr + 
    ( sum(inc_day_all_gii4[14]) +
    sum(inc_day_all_gii[14]) +
    sum(inc_day_all_gi3[14]) +
    sum(inc_day_all_gi[14]) ) * hosp_rate_elder  
    
    
  
  update(new_cases_week_gi)<- new_cases_week_gi +  sum(inc_day_gi)
  
  update(new_cases_week_gi3)<- new_cases_week_gi3 +  sum(inc_day_gi3)
  
  update(new_cases_week_gii)<- new_cases_week_gii +  sum(inc_day_gii)
  
  update(new_cases_week_gii4)<- new_cases_week_gii4 +  sum(inc_day_gii4)
  
  
  
  update(seasonality) <-season_func
  
  ## Incidence by Year and strain
  ###############
  
  
  update(inc_year_gi3[]) <- 
    if(i==1) 
      inc_year_gi3[i]+
    (sum(n_E_I[1,1])  +
       sum(n_E2_I2[1,1])  +
       sum(n_E3_I3[1,1]) +
       sum(n_E4_I4[1,1])  +
       sum(n_E23_I23[1,1])  +
       sum(n_E24_I24[1,1])  +
       sum(n_E34_I34[1,1])  +
       sum(n_E4rd_I4rd[1,1])) else if(i==2)  
         inc_year_gi3[i]+
    sum(n_E_I[2:4,1])  +
    sum(n_E2_I2[2:4,1])  +
    sum(n_E3_I3[2:4,1]) +
    sum(n_E4_I4[2:4,1])  +
    sum(n_E23_I23[2:4,1])  +
    sum(n_E24_I24[2:4,1])  +
    sum(n_E34_I34[2:4,1])  +
    sum(n_E4rd_I4rd[2:4,1]) else if(i==3)  
      inc_year_gi3[i]+
    sum(n_E_I[5:8,1])  +
    sum(n_E2_I2[5:8,1])  +
    sum(n_E3_I3[5:8,1]) +
    sum(n_E4_I4[5:8,1])  +
    sum(n_E23_I23[5:8,1])  +
    sum(n_E24_I24[5:8,1])  +
    sum(n_E34_I34[5:8,1])  +
    sum(n_E4rd_I4rd[5:8,1]) else if (i==4) 
      inc_year_gi3[i]+
    sum(n_E_I[9:13,1])  +
    sum(n_E2_I2[9:13,1])  +
    sum(n_E3_I3[9:13,1]) +
    sum(n_E4_I4[9:13,1])  +
    sum(n_E23_I23[9:13,1])  +
    sum(n_E24_I24[9:13,1])  +
    sum(n_E34_I34[9:13,1])  +
    sum(n_E4rd_I4rd[9:13,1]) else
      inc_year_gi3[i]+
    sum(n_E_I[14,1])  +
    sum(n_E2_I2[14,1])  +
    sum(n_E3_I3[14,1]) +
    sum(n_E4_I4[14,1])  +
    sum(n_E23_I23[14,1])  +
    sum(n_E24_I24[14,1])  +
    sum(n_E34_I34[14,1])  +
    sum(n_E4rd_I4rd[14,1])
  
  
  
  
  
  update(inc_year_gi[]) <-  
    if(i==1) 
      inc_year_gi[i]+
    sum(n_E_I[1,2] ) +
    sum(n_E1_I1[1,1])  +
    sum(n_E3_I3[1,2]) +
    sum(n_E4_I4[1,2])  +
    sum(n_E13_I13[1,1])  +
    sum(n_E14_I14[1,1]) +
    sum(n_E34_I34[1,2])  +
    sum(n_E4rd_I4rd[1,2])  else if(i==2)  
      inc_year_gi[i]+
    sum(n_E_I[2:4,2] ) +
    sum(n_E1_I1[2:4,1])  +
    sum(n_E3_I3[2:4,2]) +
    sum(n_E4_I4[2:4,2])  +
    sum(n_E13_I13[2:4,1])  +
    sum(n_E14_I14[2:4,1]) +
    sum(n_E34_I34[2:4,2])  +
    sum(n_E4rd_I4rd[2:4,2]) else if(i==3) 
      inc_year_gi[i]+
    sum(n_E_I[5:8,2] ) +
    sum(n_E1_I1[5:8,1])  +
    sum(n_E3_I3[5:8,2]) +
    sum(n_E4_I4[5:8,2])  +
    sum(n_E13_I13[5:8,1])  +
    sum(n_E14_I14[5:8,1]) +
    sum(n_E34_I34[5:8,2])  +
    sum(n_E4rd_I4rd[5:8,2]) else if (i==4) 
      inc_year_gi[i]+
    sum(n_E_I[9:13,2] ) +
    sum(n_E1_I1[9:13,1])  +
    sum(n_E3_I3[9:13,2]) +
    sum(n_E4_I4[9:13,2])  +
    sum(n_E13_I13[9:13,1])  +
    sum(n_E14_I14[9:13,1]) +
    sum(n_E34_I34[9:13,2])  +
    sum(n_E4rd_I4rd[9:13,2]) else
      inc_year_gi[i]+
    sum(n_E_I[14,2] ) +
    sum(n_E1_I1[14,1])  +
    sum(n_E3_I3[14,2]) +
    sum(n_E4_I4[14,2])  +
    sum(n_E13_I13[14,1])  +
    sum(n_E14_I14[14,1]) +
    sum(n_E34_I34[14,2])  +
    sum(n_E4rd_I4rd[14,2]) 
  
  
  
  update(inc_year_gii4[]) <- 
    if(i==1) 
      inc_year_gii4[i]+
    sum(n_E_I[1,3])  +
    sum(n_E1_I1[1,2])  +
    sum(n_E2_I2[1,2])  +
    sum(n_E4_I4[1,3])  +
    sum(n_E12_I12[1,1])  +
    sum(n_E14_I14[1,2]) +
    sum(n_E24_I24[1,2])  +
    sum(n_E4rd_I4rd[1,3]) else if(i==2)  
      inc_year_gii4[i]+
    sum(n_E_I[2:4,3])  +
    sum(n_E1_I1[2:4,2])  +
    sum(n_E2_I2[2:4,2])  +
    sum(n_E4_I4[2:4,3])  +
    sum(n_E12_I12[2:4,1])  +
    sum(n_E14_I14[2:4,2]) +
    sum(n_E24_I24[2:4,2])  +
    sum(n_E4rd_I4rd[2:4,3]) else if(i==3)  
      inc_year_gii4[i]+
    sum(n_E_I[5:8,3])  +
    sum(n_E1_I1[5:8,2])  +
    sum(n_E2_I2[5:8,2])  +
    sum(n_E4_I4[5:8,3])  +
    sum(n_E12_I12[5:8,1])  +
    sum(n_E14_I14[5:8,2]) +
    sum(n_E24_I24[5:8,2])  +
    sum(n_E4rd_I4rd[5:8,3]) else if (i==4) 
      inc_year_gii4[i]+
    sum(n_E_I[9:13,3])  +
    sum(n_E1_I1[9:13,2])  +
    sum(n_E2_I2[9:13,2])  +
    sum(n_E4_I4[9:13,3])  +
    sum(n_E12_I12[9:13,1])  +
    sum(n_E14_I14[9:13,2]) +
    sum(n_E24_I24[9:13,2])  +
    sum(n_E4rd_I4rd[9:13,3]) else
      inc_year_gii4[i]+
    sum(n_E_I[14,3])  +
    sum(n_E1_I1[14,2])  +
    sum(n_E2_I2[14,2])  +
    sum(n_E4_I4[14,3])  +
    sum(n_E12_I12[14,1])  +
    sum(n_E14_I14[14,2]) +
    sum(n_E24_I24[14,2])  +
    sum(n_E4rd_I4rd[14,3])
  
  update(inc_year_gii[]) <-  
    if(i==1) 
      inc_year_gii[i]+
    sum(n_E_I[1,4])  +
    sum(n_E1_I1[1,3])  +
    sum(n_E2_I2[1,3])  +
    sum(n_E3_I3[1,3]) +
    sum(n_E12_I12[1,2])  +
    sum(n_E13_I13[1,2])  +
    sum(n_E23_I23[1,2])  +
    sum(n_E4rd_I4rd[1,4]) else if(i==2)  
      inc_year_gii[i]+
    sum(n_E_I[2:4,4])  +
    sum(n_E1_I1[2:4,3])  +
    sum(n_E2_I2[2:4,3])  +
    sum(n_E3_I3[2:4,3]) +
    sum(n_E12_I12[2:4,2])  +
    sum(n_E13_I13[2:4,2])  +
    sum(n_E23_I23[2:4,2])  +
    sum(n_E4rd_I4rd[2:4,4]) else if(i==3)  
      inc_year_gii[i]+
    sum(n_E_I[5:8,4])  +
    sum(n_E1_I1[5:8,3])  +
    sum(n_E2_I2[5:8,3])  +
    sum(n_E3_I3[5:8,3]) +
    sum(n_E12_I12[5:8,2])  +
    sum(n_E13_I13[5:8,2])  +
    sum(n_E23_I23[5:8,2])  +
    sum(n_E4rd_I4rd[5:8,4]) else if (i==4) 
      inc_year_gii[i]+
    sum(n_E_I[9:13,4])  +
    sum(n_E1_I1[9:13,3])  +
    sum(n_E2_I2[9:13,3])  +
    sum(n_E3_I3[9:13,3]) +
    sum(n_E12_I12[9:13,2])  +
    sum(n_E13_I13[9:13,2])  +
    sum(n_E23_I23[9:13,2])  +
    sum(n_E4rd_I4rd[9:13,4])  else
      inc_year_gii[i]+
    sum(n_E_I[14,4])  +
    sum(n_E1_I1[14,3])  +
    sum(n_E2_I2[14,3])  +
    sum(n_E3_I3[14,3]) +
    sum(n_E12_I12[14,2])  +
    sum(n_E13_I13[14,2])  +
    sum(n_E23_I23[14,2])  +
    sum(n_E4rd_I4rd[14,4])
  
  
  ## compute at risk population ( match IDD2 data)
  update(pop_by4age[]) <-
    (if(i==1)  sum(N_byage[1]) else
      if(i==2)  sum(N_byage[2:4])  else 
        if(i==3)  sum(N_byage[5:8])  else 
          if (i==4)  sum(N_byage[9:13])  else
            sum(N_byage[14]) ) 
  
  ## compute prevalence og GII4
  
  update(seroprev_gii4[]) <- prev_byage_gii4[i]/N_byage[i]
  
  
  update(sus_gi3) <- 
    sum(S)+
    sum(R[,3:4])+
    sum(R3[,3])+
    sum(R4[,3])
  
  update(sus_cross_gi3) <- 
    sum(R[,2])+
    sum(R2[,2])+
    sum(R2[,3])+
    sum(R3[,2])+
    sum(R4[,2])+
    sum(R23[,2])+
    sum(R24[,2])+
    sum(R34[,2])
  
  update(sus_re_gi3) <- 
    sum(R[,1])+
    sum(R2[,1])+
    sum(R3[,1])+
    sum(R4[,1])+
    sum(R23[,1])+
    sum(R24[,1])+
    sum(R34[,1])+
    sum(R4rd[,1])
  
  
  update(sus_gi) <- 
    sum(S)+
    sum(R[,3:4])+
    sum(R3[,3])+
    sum(R4[,3])
  
  update(sus_cross_gi) <- 
    sum(R[,1])+
    sum(R1[,2])+
    sum(R1[,3])+
    sum(R3[,1])+
    sum(R4[,1])+
    sum(R13[,2])+
    sum(R14[,2])+
    sum(R34[,1])
  
  
  update(sus_re_gi)  <- 
    sum(R[,2])+
    sum(R1[,1])+
    sum(R3[,2])+
    sum(R4[,2])+
    sum(R13[,1])+
    sum(R14[,1])+
    sum(R34[,2])+ #################################################################
  sum(R4rd[,2])
  
  
  update(sus_gii4) <- 
    sum(S)+
    sum(R[,1:2])+
    sum(R1[,1])+
    sum(R2[,1])
  
  update(sus_cross_gii4) <-   
    sum(R[,4])+
    sum(R1[,3])+
    sum(R2[,3])+
    sum(R4[,1:2])+
    sum(R12[,2])+
    sum(R14[,1])+
    sum(R24[,1])
  
  update(sus_re_gii4) <- 
    sum(R[,3])+
    sum(R1[,2])+
    sum(R2[,2])+
    sum(R4[,3])+
    sum(R12[,1])+
    sum(R14[,2])+
    sum(R24[,2])+
    sum(R4rd[,3])
  
  
  update(sus_gii) <- 
    sum(S)+
    sum(R[,1:2])+
    sum(R1[,1])+
    sum(R2[,1])
  
  update(sus_cross_gii) <- 
    sum(R[,3])+
    sum(R1[,2])+
    sum(R2[,2])+
    sum(R3[,1:2])+
    sum(R12[,1])+
    sum(R13[,1])+
    sum(R23[,1])
  
  update(sus_re_gii) <-
    sum(R[,4])+
    sum(R1[,3])+
    sum(R2[,3])+
    sum(R3[,3])+
    sum(R12[,2])+
    sum(R13[,2])+
    sum(R23[,2])+
    sum(R4rd[,4])
  
  
  
  update(covid_timeline)<-contact_matrix[5,5]
  
  
  
  # Compare built-in function -----------------------------------------------
  
  ## load data
  cases_a1<-data()
  cases_a2<-data()
  cases_a3<-data()
  cases_a4<-data()
  cases_a5<-data()
  a0_event<-data()
  a5_event<-data()
  a15_event<-data()
  a65_event<-data()
  reported<-data()
  reported_gi3 <-data()
  reported_gi  <-data()
  reported_gii4<-data()
  reported_gii <-data()
  sero1<-data()
  sero2<-data()
  sero3<-data()
  sero4<-data()
  sero5<-data()
  sero6<-data()
  
  a0_size<-56308/12
  a5_size<-56308/12
  a15_size<-56308/12
  a65_size<-56308/12
  
  sero1_n<-103
  sero2_n<-107
  sero3_n<-121
  sero4_n<-124
  sero5_n<-122
  sero6_n<-109
  
  noise <- Exponential(1e+06)
  
  ## Reported by age outcomes
  rep_all<-sum(infections_day_gi3)+
    sum(infections_day_gi)+
    sum(infections_day_gii4)+
    sum(infections_day_gii)
  
  
  rep_a0<- (
    infections_day_gi3[1]+
      infections_day_gi[1]+
      infections_day_gii4[1]+
      infections_day_gii[1])/rep_all
  
  rep_a5<- 
    (infections_day_gi3[2]+
       infections_day_gi[2]+
       infections_day_gii4[2]+
       infections_day_gii[2])/rep_all
  
  rep_a15<- 
    (infections_day_gi3[3]+
       infections_day_gi[3]+
       infections_day_gii4[3]+
       infections_day_gii[3])/rep_all
  
  rep_a65<- 
    (infections_day_gi3[4]+
       infections_day_gi[4]+
       infections_day_gii4[4]+
       infections_day_gii[4])/rep_all
  
  # IID2 incidence -------------------------------------------------------
  cases_a1 ~ Poisson(1000 * (inc_year_gi3[1]+inc_year_gi[1]+inc_year_gii4[1]+inc_year_gii[1])/pop_by4age[1] + noise)
  cases_a2 ~ Poisson(1000 * (inc_year_gi3[2]+inc_year_gi[2]+inc_year_gii4[2]+inc_year_gii[2])/pop_by4age[2] + noise)
  cases_a3 ~ Poisson(1000 * (inc_year_gi3[3]+inc_year_gi[3]+inc_year_gii4[3]+inc_year_gii[3])/pop_by4age[3] + noise)
  cases_a4 ~ Poisson(1000 * (inc_year_gi3[4]+inc_year_gi[4]+inc_year_gii4[4]+inc_year_gii[4])/pop_by4age[4] + noise)
  cases_a5 ~ Poisson(1000 * (inc_year_gi3[5]+inc_year_gi[5]+inc_year_gii4[5]+inc_year_gii[5])/pop_by4age[5] + noise)
  
  # SGSS weekly reported series by strain-------------------------------------------------------
  reported ~ Poisson(sum(reported_wk) + noise) 
  
  
  # SGSS weekly reported series by strain-------------------------------------------------------
  reported_gi3  ~ Poisson(reported_wk_gi3 + noise)#/52
  reported_gi   ~ Poisson(reported_wk_gi  + noise)#/52
  reported_gii4 ~ Poisson(reported_wk_gii4+ noise)#/52
  reported_gii  ~ Poisson(reported_wk_gii + noise)#/52
  
  
  # Seroprev GII4 children 1 to 7 -------------------------------------------------------
  
  sero1 ~ Binomial(sero1_n, seroprev_gii4[2] + noise)
  sero2 ~ Binomial(sero2_n, seroprev_gii4[3] + noise)
  sero3 ~ Binomial(sero3_n, seroprev_gii4[4] + noise)
  sero4 ~ Binomial(sero4_n, seroprev_gii4[5] + noise)
  sero5 ~ Binomial(sero5_n, seroprev_gii4[6] + noise)
  sero6 ~ Binomial(sero6_n, seroprev_gii4[7] + noise)
  
  
  # SGSS age reporting fractions -------------------------------------------------------
  
  a0_event  ~ Binomial(a0_size, rep_a0 + noise)
  a5_event  ~ Binomial(a5_size, rep_a5 + noise)
  a15_event ~ Binomial(a15_size, rep_a15 + noise)
  a65_event ~ Binomial(a65_size, rep_a65 + noise)
  
  # Individual probabilities of transition: -------------------------------------------------------
  
  p_mu[]     <- 1 - exp(-mu[i] * dt) # mortality
  p_cfr[]    <- 1 - exp(-cfr[i] * dt) # mortality
  p_aging[]  <- 1 - exp(-aging_vec[i] * dt)
  p_MS[]     <- 1 - exp(-(1/maternalAB) * dt)  # M to S
  p_SE[,]    <- 1 - exp(- lambda[i,j] * dt) # S to E all
  p_EI       <- 1 - exp(-(1/epsilon) * dt) # E to I
  p_IA_5     <- 1 - exp(-(1/theta_5) * dt) # I to A
  p_IA_5p    <- 1 - exp(-(1/theta_5p) * dt) # I to A
  p_AR       <- 1 - exp(-(1/sigma) * dt) # A to R
  p_wane       <- 1 - exp(- (1/(imm_yr*365)) * dt) # R to S
  p_RS2      <- 1 - exp(- (1/(imm_yr*imm_fac*365)) * dt) # R to S
  p_vacc[]   <- 1 - exp(- vaccination_coverage[i] * dt) # mortality
  
  N_byage[]<-(   
      G[i] +
    M[i] +  
      S[i] + 
      sum(E[i,])+
      sum(I[i,])+
      sum(A[i,])+
      sum(R[i,])+
      sum(E1[i,])+
      sum(I1[i,])+
      sum(A1[i,])+
      sum(R1[i,])+
      sum(E2[i,])+
      sum(I2[i,])+
      sum(A2[i,])+
      sum(R2[i,])+
      sum(E3[i,])+
      sum(I3[i,])+
      sum(A3[i,])+
      sum(R3[i,])+
      sum(E4[i,])+
      sum(I4[i,])+
      sum(A4[i,])+
      sum(R4[i,])+
      sum(E12[i,])+
      sum(I12[i,])+
      sum(A12[i,])+
      sum(R12[i,])+
      sum(E13[i,])+
      sum(I13[i,])+
      sum(A13[i,])+
      sum(R13[i,])+
      sum(E14[i,])+
      sum(I14[i,])+
      sum(A14[i,])+
      sum(R14[i,])+
      sum(E23[i,])+
      sum(I23[i,])+
      sum(A23[i,])+
      sum(R23[i,])+
      sum(E24[i,])+
      sum(I24[i,])+
      sum(A24[i,])+
      sum(R24[i,])+
      sum(E34[i,])+
      sum(I34[i,])+
      sum(A34[i,])+
      sum(R34[i,])+
      sum(E4rd[i,])+
      sum(I4rd[i,])+
      sum(A4rd[i,])+
      sum(R4rd[i,]))
  
  
  
  prev_byage_gii4[]<-(M[i] +
                        I[i,3]+
                        A[i,3]+
                        R[i,3]+
                        I1[i,2]+
                        A1[i,2]+
                        R1[i,2]+
                        I2[i,2]+
                        A2[i,2]+
                        R2[i,2]+
                        sum(I3[i,])+
                        sum(A3[i,])+
                        sum(R3[i,])+
                        I4[i,3]+
                        A4[i,3]+
                        R4[i,3]+
                        I12[i,1]+
                        A12[i,1]+
                        R12[i,1]+
                        sum(I13[i,])+
                        sum(A13[i,])+
                        sum(R13[i,])+
                        I14[i,2]+
                        A14[i,2]+
                        R14[i,2]+
                        sum(I23[i,])+
                        sum(A23[i,])+
                        sum(R23[i,])+
                        I24[i,2]+
                        A24[i,2]+
                        R24[i,2]+
                        sum(I34[i,])+
                        sum(A34[i,])+
                        sum(R34[i,])+
                        I4rd[i,3]+
                        A4rd[i,3]+
                        R4rd[i,3])
  
  prev_byage_all[]<-(
    M[i]+
      sum(I[i,])+
      sum(A[i,])+
      sum(R[i,])+
      sum(I1[i,])+
      sum(A1[i,])+
      sum(R1[i,])+
      sum(I2[i,])+
      sum(A2[i,])+
      sum(R2[i,])+
      sum(I3[i,])+
      sum(A3[i,])+
      sum(R3[i,])+
      sum(I4[i,])+
      sum(A4[i,])+
      sum(R4[i,])+
      sum(I12[i,])+
      sum(A12[i,])+
      sum(R12[i,])+
      sum(I13[i,])+
      sum(A13[i,])+
      sum(R13[i,])+
      sum(I14[i,])+
      sum(A14[i,])+
      sum(R14[i,])+
      sum(I23[i,])+
      sum(A23[i,])+
      sum(R23[i,])+
      sum(I24[i,])+
      sum(A24[i,])+
      sum(R24[i,])+
      sum(I34[i,])+
      sum(A34[i,])+
      sum(R34[i,])+
      sum(I4rd[i,])+
      sum(A4rd[i,])+
      sum(R4rd[i,])
  )
  
  
  
  
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
  c1_ij[] <- (I[i,1]+I2[i,1] + 
                I3[i,1]+I4[i,1] +
                I23[i,1]+I24[i,1] + 
                I34[i,1]+I4rd[i,1] + 
                rr_inf_asymp * (A[i,1]+A2[i,1] + 
                                  A3[i,1]+A4[i,1] +
                                  A23[i,1]+A24[i,1] + 
                                  A34[i,1]+A4rd[i,1] +
                                  E[i,1]+E2[i,1] + 
                                  E3[i,1]+E4[i,1] +
                                  E23[i,1]+E24[i,1] + 
                                  E34[i,1]+E4rd[i,1]) ) * age_RR[i] 
  
  
  
  # infectives GI
  c2_ij[] <- (I[i,2]+
                I1[i,1]+I3[i,2]+I4[i,2]+
                I13[i,1]+I14[i,1]+I34[i,2]+I4rd[i,2] +
                rr_inf_asymp * (A[i,2]+
                                  A1[i,1]+A3[i,2]+A4[i,2]+
                                  A13[i,1]+A14[i,1]+ A34[i,2]+A4rd[i,2]+
                                  E[i,2]+
                                  E1[i,1]+E3[i,2]+E4[i,2]+
                                  E13[i,1]+E14[i,1]+ E34[i,2]+
                                  E4rd[i,2]))   * age_RR[i] 
  
  
  # infectives GII4
  c3_ij[] <- (I[i,3]+
                I1[i,2]+I2[i,2]+I4[i,3]+
                I12[i,1]+I14[i,2]+I24[i,2]+
                I4rd[i,3]+ 
                rr_inf_asymp *(A[i,3]+
                                 A1[i,2]+A2[i,2]+A4[i,3]+
                                 A12[i,1]+A14[i,2]+A24[i,2]+
                                 A4rd[i,3]+
                                 E[i,3]+
                                 E1[i,2]+E2[i,2]+E4[i,3]+
                                 E12[i,1]+E14[i,2]+ E24[i,2]+
                                 E4rd[i,3]) )  * age_RR[i] 
  
  
  
  
  c4_ij[] <- (I[i,4]+
                I1[i,3]+I2[i,3]+I3[i,3]+
                I12[i,2]+I13[i,2]+I23[i,2]+
                I4rd[i,4]+ 
                rr_inf_asymp *(A[i,4]+
                                 A1[i,3]+A2[i,3]+A3[i,3]+
                                 A12[i,2]+A13[i,2]+A23[i,2]+
                                 A4rd[i,4]+
                                 E[i,4]+
                                 E1[i,3]+E2[i,3]+E3[i,3]+
                                 E12[i,2]+E13[i,2]+ E23[i,2]+
                                 E4rd[i,4]) )   * age_RR[i] 
  
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
  
  
  update(foi_gi3[]) <-rel_foi_strain[i,1]
  update(foi_gi[])  <-rel_foi_strain[i,2]
  update(foi_gii4[])<-rel_foi_strain[i,3]
  update(foi_gii[]) <-rel_foi_strain[i,4]
  
  
  
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
  n_ageoE1[,] <- Binomial(E1[i,j] ,  p_aging[i])
  n_ageoI1[,] <- Binomial(I1[i,j] ,  p_aging[i])
  n_ageoA1[,] <- Binomial(A1[i,j] ,  p_aging[i])
  n_ageoR1[,] <- Binomial(R1[i,j] ,  p_aging[i])
  n_ageoE2[,] <- Binomial(E2[i,j] ,  p_aging[i])
  n_ageoI2[,] <- Binomial(I2[i,j] ,  p_aging[i])
  n_ageoA2[,] <- Binomial(A2[i,j] ,  p_aging[i])
  n_ageoR2[,] <- Binomial(R2[i,j] ,  p_aging[i])
  n_ageoE3[,] <- Binomial(E3[i,j] ,  p_aging[i])
  n_ageoI3[,] <- Binomial(I3[i,j] ,  p_aging[i])
  n_ageoA3[,] <- Binomial(A3[i,j] ,  p_aging[i])
  n_ageoR3[,] <- Binomial(R3[i,j] ,  p_aging[i])
  n_ageoE4[,] <- Binomial(E4[i,j] ,  p_aging[i])
  n_ageoI4[,] <- Binomial(I4[i,j] ,  p_aging[i])
  n_ageoA4[,] <- Binomial(A4[i,j] ,  p_aging[i])
  n_ageoR4[,] <- Binomial(R4[i,j] ,  p_aging[i])
  n_ageoE12[,] <- Binomial(E12[i,j] ,  p_aging[i])
  n_ageoI12[,] <- Binomial(I12[i,j] ,  p_aging[i])
  n_ageoA12[,] <- Binomial(A12[i,j] ,  p_aging[i])
  n_ageoR12[,] <- Binomial(R12[i,j] ,  p_aging[i])
  n_ageoE13[,] <- Binomial(E13[i,j] ,  p_aging[i])
  n_ageoI13[,] <- Binomial(I13[i,j] ,  p_aging[i])
  n_ageoA13[,] <- Binomial(A13[i,j] ,  p_aging[i])
  n_ageoR13[,] <- Binomial(R13[i,j] ,  p_aging[i])
  n_ageoE14[,] <- Binomial(E14[i,j] ,  p_aging[i])
  n_ageoI14[,] <- Binomial(I14[i,j] ,  p_aging[i])
  n_ageoA14[,] <- Binomial(A14[i,j] ,  p_aging[i])
  n_ageoR14[,] <- Binomial(R14[i,j] ,  p_aging[i])
  n_ageoE23[,] <- Binomial(E23[i,j] ,  p_aging[i])
  n_ageoI23[,] <- Binomial(I23[i,j] ,  p_aging[i])
  n_ageoA23[,] <- Binomial(A23[i,j] ,  p_aging[i])
  n_ageoR23[,] <- Binomial(R23[i,j] ,  p_aging[i])
  n_ageoE24[,] <- Binomial(E24[i,j] ,  p_aging[i])
  n_ageoI24[,] <- Binomial(I24[i,j] ,  p_aging[i])
  n_ageoA24[,] <- Binomial(A24[i,j] ,  p_aging[i])
  n_ageoR24[,] <- Binomial(R24[i,j] ,  p_aging[i])
  n_ageoE34[,] <- Binomial(E34[i,j] ,  p_aging[i])
  n_ageoI34[,] <- Binomial(I34[i,j] ,  p_aging[i])
  n_ageoA34[,] <- Binomial(A34[i,j] ,  p_aging[i])
  n_ageoR34[,] <- Binomial(R34[i,j] ,  p_aging[i])
  n_ageoE4rd[,] <- Binomial(E4rd[i,j] ,  p_aging[i])
  n_ageoI4rd[,] <- Binomial(I4rd[i,j] ,  p_aging[i])
  n_ageoA4rd[,] <- Binomial(A4rd[i,j] ,  p_aging[i])
  n_ageoR4rd[,] <- Binomial(R4rd[i,j] ,  p_aging[i])
  
  
  
  
  #Age in
  n_ageiG[] <- if (i>1) Binomial(G[i-1] ,p_aging  [i-1]) else 0
  n_ageiM[] <- if (i>1) Binomial(M[i-1] ,p_aging  [i-1]) else 0
  n_ageiS[] <- if (i>1) Binomial(S[i-1] ,p_aging  [i-1]) else 0
  n_ageiE[,] <- if (i>1) Binomial(E[i-1,j] , p_aging[i-1])else 0
  n_ageiI[,] <- if (i>1) Binomial(I[i-1,j] , p_aging[i-1])else 0
  n_ageiA[,] <- if (i>1) Binomial(A[i-1,j] , p_aging[i-1])else 0
  n_ageiR[,] <- if (i>1) Binomial(R[i-1,j] , p_aging[i-1])else 0
  n_ageiE1[,] <- if (i>1) Binomial(E1[i-1,j] , p_aging[i-1])else 0
  n_ageiI1[,] <- if (i>1) Binomial(I1[i-1,j] , p_aging[i-1])else 0
  n_ageiA1[,] <- if (i>1) Binomial(A1[i-1,j] , p_aging[i-1])else 0
  n_ageiR1[,] <- if (i>1) Binomial(R1[i-1,j] , p_aging[i-1])else 0
  n_ageiE2[,] <- if (i>1) Binomial(E2[i-1,j] , p_aging[i-1])else 0
  n_ageiI2[,] <- if (i>1) Binomial(I2[i-1,j] , p_aging[i-1])else 0
  n_ageiA2[,] <- if (i>1) Binomial(A2[i-1,j] , p_aging[i-1])else 0
  n_ageiR2[,] <- if (i>1) Binomial(R2[i-1,j] , p_aging[i-1])else 0
  n_ageiE3[,] <- if (i>1) Binomial(E3[i-1,j] , p_aging[i-1])else 0
  n_ageiI3[,] <- if (i>1) Binomial(I3[i-1,j] , p_aging[i-1])else 0
  n_ageiA3[,] <- if (i>1) Binomial(A3[i-1,j] , p_aging[i-1])else 0
  n_ageiR3[,] <- if (i>1) Binomial(R3[i-1,j] , p_aging[i-1])else 0
  n_ageiE4[,] <- if (i>1) Binomial(E4[i-1,j] , p_aging[i-1])else 0
  n_ageiI4[,] <- if (i>1) Binomial(I4[i-1,j] , p_aging[i-1])else 0
  n_ageiA4[,] <- if (i>1) Binomial(A4[i-1,j] , p_aging[i-1])else 0
  n_ageiR4[,] <- if (i>1) Binomial(R4[i-1,j] , p_aging[i-1])else 0
  n_ageiE12[,] <- if (i>1) Binomial(E12[i-1,j] , p_aging[i-1])else 0
  n_ageiI12[,] <- if (i>1) Binomial(I12[i-1,j] , p_aging[i-1])else 0
  n_ageiA12[,] <- if (i>1) Binomial(A12[i-1,j] , p_aging[i-1])else 0
  n_ageiR12[,] <- if (i>1) Binomial(R12[i-1,j] , p_aging[i-1])else 0
  n_ageiE13[,] <- if (i>1) Binomial(E13[i-1,j] , p_aging[i-1])else 0
  n_ageiI13[,] <- if (i>1) Binomial(I13[i-1,j] , p_aging[i-1])else 0
  n_ageiA13[,] <- if (i>1) Binomial(A13[i-1,j] , p_aging[i-1])else 0
  n_ageiR13[,] <- if (i>1) Binomial(R13[i-1,j] , p_aging[i-1])else 0
  n_ageiE14[,] <- if (i>1) Binomial(E14[i-1,j] , p_aging[i-1])else 0
  n_ageiI14[,] <- if (i>1) Binomial(I14[i-1,j] , p_aging[i-1])else 0
  n_ageiA14[,] <- if (i>1) Binomial(A14[i-1,j] , p_aging[i-1])else 0
  n_ageiR14[,] <- if (i>1) Binomial(R14[i-1,j] , p_aging[i-1])else 0
  n_ageiE23[,] <- if (i>1) Binomial(E23[i-1,j] , p_aging[i-1])else 0
  n_ageiI23[,] <- if (i>1) Binomial(I23[i-1,j] , p_aging[i-1])else 0
  n_ageiA23[,] <- if (i>1) Binomial(A23[i-1,j] , p_aging[i-1])else 0
  n_ageiR23[,] <- if (i>1) Binomial(R23[i-1,j] , p_aging[i-1])else 0
  n_ageiE24[,] <- if (i>1) Binomial(E24[i-1,j] , p_aging[i-1])else 0
  n_ageiI24[,] <- if (i>1) Binomial(I24[i-1,j] , p_aging[i-1])else 0
  n_ageiA24[,] <- if (i>1) Binomial(A24[i-1,j] , p_aging[i-1])else 0
  n_ageiR24[,] <- if (i>1) Binomial(R24[i-1,j] , p_aging[i-1])else 0
  n_ageiE34[,] <- if (i>1) Binomial(E34[i-1,j] , p_aging[i-1])else 0
  n_ageiI34[,] <- if (i>1) Binomial(I34[i-1,j] , p_aging[i-1])else 0
  n_ageiA34[,] <- if (i>1) Binomial(A34[i-1,j] , p_aging[i-1])else 0
  n_ageiR34[,] <- if (i>1) Binomial(R34[i-1,j] , p_aging[i-1])else 0
  n_ageiE4rd[,] <- if (i>1) Binomial(E4rd[i-1,j] , p_aging[i-1])else 0
  n_ageiI4rd[,] <- if (i>1) Binomial(I4rd[i-1,j] , p_aging[i-1])else 0
  n_ageiA4rd[,] <- if (i>1) Binomial(A4rd[i-1,j] , p_aging[i-1])else 0
  n_ageiR4rd[,] <- if (i>1) Binomial(R4rd[i-1,j] , p_aging[i-1])else 0
  
  
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
  n_muE1[, ] <- Binomial(E1[i, j] - n_ageoE1[i, j], p_mu[i])
  n_muI1[, ] <- Binomial(I1[i, j] - n_ageoI1[i, j], p_mu[i])
  n_muA1[, ] <- Binomial(A1[i, j] - n_ageoA1[i, j], p_mu[i])
  n_muR1[, ] <- Binomial(R1[i, j] - n_ageoR1[i, j], p_mu[i])
  n_muE2[, ] <- Binomial(E2[i, j] - n_ageoE2[i, j], p_mu[i])
  n_muI2[, ] <- Binomial(I2[i, j] - n_ageoI2[i, j], p_mu[i])
  n_muA2[, ] <- Binomial(A2[i, j] - n_ageoA2[i, j], p_mu[i])
  n_muR2[, ] <- Binomial(R2[i, j] - n_ageoR2[i, j], p_mu[i])
  n_muE3[, ] <- Binomial(E3[i, j] - n_ageoE3[i, j], p_mu[i])
  n_muI3[, ] <- Binomial(I3[i, j] - n_ageoI3[i, j], p_mu[i])
  n_muA3[, ] <- Binomial(A3[i, j] - n_ageoA3[i, j], p_mu[i])
  n_muR3[, ] <- Binomial(R3[i, j] - n_ageoR3[i, j], p_mu[i])
  n_muE4[, ] <- Binomial(E4[i, j] - n_ageoE4[i, j], p_mu[i])
  n_muI4[, ] <- Binomial(I4[i, j] - n_ageoI4[i, j], p_mu[i])
  n_muA4[, ] <- Binomial(A4[i, j] - n_ageoA4[i, j], p_mu[i])
  n_muR4[, ] <- Binomial(R4[i, j] - n_ageoR4[i, j], p_mu[i])
  n_muE12[, ] <- Binomial(E12[i, j] - n_ageoE12[i, j], p_mu[i])
  n_muI12[, ] <- Binomial(I12[i, j] - n_ageoI12[i, j], p_mu[i])
  n_muA12[, ] <- Binomial(A12[i, j] - n_ageoA12[i, j], p_mu[i])
  n_muR12[, ] <- Binomial(R12[i, j] - n_ageoR12[i, j], p_mu[i])
  n_muE13[, ] <- Binomial(E13[i, j] - n_ageoE13[i, j], p_mu[i])
  n_muI13[, ] <- Binomial(I13[i, j] - n_ageoI13[i, j], p_mu[i])
  n_muA13[, ] <- Binomial(A13[i, j] - n_ageoA13[i, j], p_mu[i])
  n_muR13[, ] <- Binomial(R13[i, j] - n_ageoR13[i, j], p_mu[i])
  n_muE14[, ] <- Binomial(E14[i, j] - n_ageoE14[i, j], p_mu[i])
  n_muI14[, ] <- Binomial(I14[i, j] - n_ageoI14[i, j], p_mu[i])
  n_muA14[, ] <- Binomial(A14[i, j] - n_ageoA14[i, j], p_mu[i])
  n_muR14[, ] <- Binomial(R14[i, j] - n_ageoR14[i, j], p_mu[i])
  n_muE23[, ] <- Binomial(E23[i, j] - n_ageoE23[i, j], p_mu[i])
  n_muI23[, ] <- Binomial(I23[i, j] - n_ageoI23[i, j], p_mu[i])
  n_muA23[, ] <- Binomial(A23[i, j] - n_ageoA23[i, j], p_mu[i])
  n_muR23[, ] <- Binomial(R23[i, j] - n_ageoR23[i, j], p_mu[i])
  n_muE24[, ] <- Binomial(E24[i, j] - n_ageoE24[i, j], p_mu[i])
  n_muI24[, ] <- Binomial(I24[i, j] - n_ageoI24[i, j], p_mu[i])
  n_muA24[, ] <- Binomial(A24[i, j] - n_ageoA24[i, j], p_mu[i])
  n_muR24[, ] <- Binomial(R24[i, j] - n_ageoR24[i, j], p_mu[i])
  n_muE34[, ] <- Binomial(E34[i, j] - n_ageoE34[i, j], p_mu[i])
  n_muI34[, ] <- Binomial(I34[i, j] - n_ageoI34[i, j], p_mu[i])
  n_muA34[, ] <- Binomial(A34[i, j] - n_ageoA34[i, j], p_mu[i])
  n_muR34[, ] <- Binomial(R34[i, j] - n_ageoR34[i, j], p_mu[i])
  n_muE4rd[, ] <- Binomial(E4rd[i, j] - n_ageoE4rd[i, j], p_mu[i])
  n_muI4rd[, ] <- Binomial(I4rd[i, j] - n_ageoI4rd[i, j], p_mu[i])
  n_muA4rd[, ] <- Binomial(A4rd[i, j] - n_ageoA4rd[i, j], p_mu[i])
  n_muR4rd[, ] <- Binomial(R4rd[i, j] - n_ageoR4rd[i, j], p_mu[i])
  
  
  # Vaccine switches --------------------------------------------------------
  
  
  switch_vacc[]  <- if (time>t_zero) vacc_switch_on[i] else 0
  switch_vacc2[] <- if (time>t_campaign_stop) 1-campaign_switch[i] else 1
  switch_vacc3[] <- switch_vacc[i]*switch_vacc2[i]
  
  dim(switch_vacc)  <- N_age
  dim(switch_vacc2) <- N_age
  dim(switch_vacc3) <- N_age
  
  
  
  # M transitions -----------------------------------------------------------
  
  
  
  n_MS[]   <- Binomial(M[i] - n_ageoM[i] - n_muM[i], p_MS[i])
  
  n_Mvacc_c[]<- Binomial(M[i] - n_ageoM[i] - n_muM[i] - n_MS[i], p_vacc[i]* switch_vacc3[i])
  
  n_Mvacc[]<- n_Mvacc_c[i] * vacc_trans
  
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
  dim(n_Mvacc_c)   <- N_age
  
  # S to E transitions ------------------------------------------------------
  
  
  n_Svacc_c[] <- Binomial(S[i] - n_ageoS[i] - n_muS[i], p_vacc[i]* switch_vacc3[i])
  
  n_Svacc[]   <- n_Svacc_c[i]*vacc_trans
  
  n_S_E_tot[] <- Binomial(S[i] - n_ageoS[i] - n_muS[i]-n_Svacc[i], sum(p_SE[i,]))
  
  
  # Multinomial draw (via nested binomial) for multiple strain infection
  n_S_E[, 1] <- Binomial(n_S_E_tot[i] , rel_foi_strain[i, 1])
  n_S_E[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial((n_S_E_tot[i] - sum(n_S_E[i, 1:(j - 1)])) , rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  
  n_S_E[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    n_S_E_tot[i] - sum(n_S_E[i, 1:3]) 
  
  dim(n_S_E_tot) <- N_age
  dim(n_Svacc)   <- N_age
  dim(n_Svacc_c)   <- N_age
  
  
  
  ## Vaccinations matrix from E (only used for accounting doses)
  n_Evacc[, ]    <-Binomial(E[i, j] - n_ageoE[i, j] - n_muE[i, j]-n_E_I[i, j],  p_vacc[i]* switch_vacc3[i])
  n_E1vacc[,]   <- Binomial(E1[i, j] - n_ageoE1[i, j] - n_muE1[i, j] - n_E1_I1[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E2vacc[,]   <- Binomial(E2[i, j] - n_ageoE2[i, j] - n_muE2[i, j] - n_E2_I2[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E3vacc[,]   <- Binomial(E3[i, j] - n_ageoE3[i, j] - n_muE3[i, j] - n_E3_I3[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E4vacc[,]   <- Binomial(E4[i, j] - n_ageoE4[i, j] - n_muE4[i, j] - n_E4_I4[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E12vacc[,]  <- Binomial(E12[i, j] - n_ageoE12[i, j] - n_muE12[i, j]-n_E12_I12[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E13vacc[,]  <- Binomial(E13[i, j] - n_ageoE13[i, j] - n_muE13[i, j]-n_E13_I13[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E14vacc[,]  <- Binomial(E14[i, j] - n_ageoE14[i, j] - n_muE14[i, j]-n_E14_I14[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E23vacc[,]  <- Binomial(E23[i, j] - n_ageoE23[i, j] - n_muE23[i, j]-n_E23_I23[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E24vacc[,]  <- Binomial(E24[i, j] - n_ageoE24[i, j] - n_muE24[i, j]-n_E24_I24[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E34vacc[,]   <- Binomial(E34[i, j] - n_ageoE34[i, j] - n_muE34[i, j]-n_E34_I34[i,j],  p_vacc[i]* switch_vacc3[i])
  n_E4rdvacc[,] <- Binomial(E4rd[i, j] - n_ageoE4rd[i, j] - n_muE4rd[i, j]-n_E4rd_I4rd[i,j],  p_vacc[i]* switch_vacc3[i])
  
  
  dim(n_Evacc	  )<-c(N_age,N_strain)
  dim(n_E1vacc	)<-c(N_age,N_strain_1)
  dim(n_E2vacc	)<-c(N_age,N_strain_1)
  dim(n_E3vacc	)<-c(N_age,N_strain_1)
  dim(n_E4vacc	)<-c(N_age,N_strain_1)
  dim(n_E12vacc	)<-c(N_age,N_strain_2)
  dim(n_E13vacc	)<-c(N_age,N_strain_2)
  dim(n_E14vacc	)<-c(N_age,N_strain_2)
  dim(n_E23vacc	)<-c(N_age,N_strain_2)
  dim(n_E24vacc	)<-c(N_age,N_strain_2)
  dim(n_E34vacc	)<-c(N_age,N_strain_2)
  dim(n_E4rdvacc)<-c(N_age,N_strain)
  
  
  ## Transition E to I
  
  
  n_E_I[, ] <- Binomial(E[i, j] - n_ageoE[i, j] - n_muE[i, j], p_EI)
  n_E_vacc[, ] <- n_Evacc[i,j] * vacc_dis 
  # n_E_vacc[, 1] <- n_Evacc[i,1] * vacc_dis 
  # n_E_vacc[, 2] <- 0
  # n_E_vacc[, 3] <- n_Evacc[i,3] * vacc_dis
  # n_E_vacc[, 4] <- 0
  
  
  
  n_E1_I1[, ]   <- Binomial(E1[i, j] - n_ageoE1[i, j] - n_muE1[i, j], p_EI)
  n_E1_vacc[,] <- n_E1vacc[i,j]* vacc_dis
  # n_E1_vacc[,1] <- 0
  # n_E1_vacc[,2] <- n_E1vacc[i,2]* vacc_dis
  # n_E1_vacc[,3] <- 0
  
  n_E2_I2[, ] <- Binomial(E2[i, j] - n_ageoE2[i, j] - n_muE2[i, j], p_EI)
  n_E2_vacc[,] <- n_E2vacc[i,j]* vacc_dis
  # n_E2_vacc[,1] <- n_E2vacc[i,1]* vacc_dis
  # n_E2_vacc[,2] <- n_E2vacc[i,2]* vacc_dis
  # n_E2_vacc[,3] <- 0
  
  n_E3_I3[, ]   <- Binomial(E3[i, j] - n_ageoE3[i, j] - n_muE3[i, j], p_EI)
  n_E3_vacc[,] <- n_E3vacc[i,j]* vacc_dis
  # n_E3_vacc[,1] <- n_E3vacc[i,1]* vacc_dis
  # n_E3_vacc[,2] <- 0
  # n_E3_vacc[,3] <- 0
  
  n_E4_I4[, ] <- Binomial(E4[i, j] - n_ageoE4[i, j] - n_muE4[i, j], p_EI)
  n_E4_vacc[,] <- n_E4vacc[i,j]* vacc_dis
  # n_E4_vacc[,1] <- n_E4vacc[i,1]* vacc_dis
  # n_E4_vacc[,2] <- 0
  # n_E4_vacc[,3] <- n_E4vacc[i,3]* vacc_dis
  
  n_E12_I12[, ]  <- Binomial(E12[i, j] - n_ageoE12[i, j] - n_muE12[i, j], p_EI)
  n_E12_vacc[,]  <- n_E12vacc[i,j]* vacc_dis
  
  # n_E12_vacc[,1] <- n_E12vacc[i,1]* vacc_dis
  # n_E12_vacc[,2] <- 0
  
  n_E13_I13[, ] <- Binomial(E13[i, j] - n_ageoE13[i, j] - n_muE13[i, j], p_EI)
  n_E13_vacc[,]  <- n_E13vacc[i,j]* vacc_dis
  # n_E13_vacc[,1] <- 0
  # n_E13_vacc[,2] <- 0
  
  n_E14_I14[, ]  <- Binomial(E14[i, j] - n_ageoE14[i, j] - n_muE14[i, j], p_EI)
  n_E14_vacc[,]  <-  n_E14vacc[i,j] * vacc_dis
  # n_E14_vacc[,1] <- 0
  # n_E14_vacc[,2] <-  n_E14vacc[i,2]* vacc_dis
  
  n_E23_I23[, ] <- Binomial(E23[i, j] - n_ageoE23[i, j] - n_muE23[i, j], p_EI)
  n_E23_vacc[,] <- n_E23vacc[i,j]* vacc_dis
  # n_E23_vacc[,1] <- n_E23vacc[i,1]* vacc_dis
  # n_E23_vacc[,2] <- 0
  
  n_E24_I24[, ]  <- Binomial(E24[i, j] - n_ageoE24[i, j] - n_muE24[i, j], p_EI)
  n_E24_vacc[,] <- n_E24vacc[i,j]* vacc_dis
  # n_E24_vacc[,1] <- n_E24vacc[i,1]* vacc_dis
  # n_E24_vacc[,2] <- n_E24vacc[i,2]* vacc_dis
  
  n_E34_I34[, ]  <- Binomial(E34[i, j] - n_ageoE34[i, j] - n_muE34[i, j], p_EI)
  n_E34_vacc[,] <- n_E34vacc[i,j]* vacc_dis
  # n_E34_vacc[,1] <- n_E34vacc[i,1]* vacc_dis
  # n_E34_vacc[,1] <- 0
  
  n_E4rd_I4rd[, ] <- Binomial(E4rd[i, j] - n_ageoE4rd[i, j] - n_muE4rd[i, j], p_EI)
  n_E4rd_vacc[,] <- n_E4rdvacc[i,j]* vacc_dis
  
  # n_E4rd_vacc[,1] <- n_E4rdvacc[i,1]* vacc_dis
  # n_E4rd_vacc[,2] <- 0
  # n_E4rd_vacc[,3] <- n_E4rdvacc[i,3]* vacc_dis
  # n_E4rd_vacc[,4] <- 0
  
  # Transitions I to A
  n_I_A[, ] <- Binomial(I[i, j] - n_ageoI[i, j] - n_muI[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I1_A1[, ] <- Binomial(I1[i, j] - n_ageoI1[i, j] - n_muI1[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I2_A2[, ] <- Binomial(I2[i, j] - n_ageoI2[i, j] - n_muI2[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I3_A3[, ] <- Binomial(I3[i, j] - n_ageoI3[i, j] - n_muI3[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I4_A4[, ] <- Binomial(I4[i, j] - n_ageoI4[i, j] - n_muI4[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I12_A12[, ] <- Binomial(I12[i, j] - n_ageoI12[i, j] - n_muI12[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I13_A13[, ] <- Binomial(I13[i, j] - n_ageoI13[i, j] - n_muI13[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I14_A14[, ] <- Binomial(I14[i, j] - n_ageoI14[i, j] - n_muI14[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I23_A23[, ] <- Binomial(I23[i, j] - n_ageoI23[i, j] - n_muI23[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I24_A24[, ] <- Binomial(I24[i, j] - n_ageoI24[i, j] - n_muI24[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I34_A34[, ] <- Binomial(I34[i, j] - n_ageoI34[i, j] - n_muI34[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  n_I4rd_A4rd[, ] <- Binomial(I4rd[i, j] - n_ageoI4rd[i, j] - n_muI4rd[i, j], if (i <= 5) p_IA_5 else p_IA_5p)
  
  
  
  
  
  # A to R transitions
  
  n_A_R[, ] <- Binomial(A[i, j] - n_ageoA[i, j] - n_muA[i, j], p_AR)
  n_A1_R1[, ] <- Binomial(A1[i, j] - n_ageoA1[i, j] - n_muA1[i, j], p_AR)
  n_A2_R2[, ] <- Binomial(A2[i, j] - n_ageoA2[i, j] - n_muA2[i, j], p_AR)
  n_A3_R3[, ] <- Binomial(A3[i, j] - n_ageoA3[i, j] - n_muA3[i, j], p_AR)
  n_A4_R4[, ] <- Binomial(A4[i, j] - n_ageoA4[i, j] - n_muA4[i, j], p_AR)
  n_A12_R12[, ] <- Binomial(A12[i, j] - n_ageoA12[i, j] - n_muA12[i, j], p_AR)
  n_A13_R13[, ] <- Binomial(A13[i, j] - n_ageoA13[i, j] - n_muA13[i, j], p_AR)
  n_A14_R14[, ] <- Binomial(A14[i, j] - n_ageoA14[i, j] - n_muA14[i, j], p_AR)
  n_A23_R23[, ] <- Binomial(A23[i, j] - n_ageoA23[i, j] - n_muA23[i, j], p_AR)
  n_A24_R24[, ] <- Binomial(A24[i, j] - n_ageoA24[i, j] - n_muA24[i, j], p_AR)
  n_A34_R34[, ] <- Binomial(A34[i, j] - n_ageoA34[i, j] - n_muA34[i, j], p_AR)
  n_A4rd_R4rd[, ] <- Binomial(A4rd[i, j] - n_ageoA4rd[i, j] - n_muA4rd[i, j], p_AR)
  
  
  
  
  ########## Transitions R 
 # create a map of cross-protection
 cross_map[,] <-if (i==1) ( 
  if(j==2) crossp_GI else 1 ) else if(i==2)(
    if(j==1) crossp_GI else 1) else if(i==3)(
      if(j==4) crossp_GII else 1) else (
        if(j==3) crossp_GII else 1)

  dim(cross_map)<-c(N_strain,N_strain)


  n_R_wane[, ] <- Binomial(R[i, j] - n_ageoR[i, j] - n_muR[i, j], p_wane)
  n_Rvacc_c[, ] <- Binomial(R[i, j] - n_ageoR[i, j] - n_muR[i, j]-n_R_wane[i,j ],p_vacc[i]* switch_vacc3[i])
  n_Rvacc[, ] <- n_Rvacc_c[i,j] * vacc_trans
  
  
  
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
  
  n_R_E1[,1] <-n_R_rel1[i,2]
  n_R_E1[,2] <-n_R_rel1[i,3]
  n_R_E1[,3] <-n_R_rel1[i,4]
  
  # R[2]
  n_R_rel2[, 1] <- Binomial(sum(n_R_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R_tot[i, 2,]) - sum(n_R_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R_tot[i, 2,]) - sum(n_R_rel2[i, 1:3])
 
  
  n_R_E2[,1] <-n_R_rel2[i,1]  
  n_R_E2[,2] <-n_R_rel2[i,3]
  n_R_E2[,3] <-n_R_rel2[i,4]

  
  # R[3]
  n_R_rel3[, 1] <- Binomial(sum(n_R_tot[i, 3,]), rel_foi_strain[i, 1])
  n_R_rel3[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R_tot[i, 3,]) - sum(n_R_rel3[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R_rel3[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R_tot[i, 3,]) - sum(n_R_rel3[i, 1:3])
  
  n_R_E3[,1] <- n_R_rel3[i,1]
  n_R_E3[,2] <- n_R_rel3[i,2]
  n_R_E3[,3] <- n_R_rel3[i,4]
  
  
  # R[4]
  n_R_rel4[, 1] <- Binomial(sum(n_R_tot[i, 4,]), rel_foi_strain[i, 1])
  n_R_rel4[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R_tot[i, 4,]) - sum(n_R_rel4[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R_rel4[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R_tot[i, 4,]) - sum(n_R_rel4[i, 1:3])
  
  n_R_E4[,1] <- n_R_rel4[i,1]
  n_R_E4[,2] <- n_R_rel4[i,2]
  n_R_E4[,3] <- n_R_rel4[i,3]
  
  dim(n_R_tot) <- c(N_age, N_strain,N_strain)
  dim(n_R_rel1) <- c(N_age, N_strain)
  dim(n_R_rel2) <- c(N_age, N_strain)
  dim(n_R_rel3) <- c(N_age, N_strain)
  dim(n_R_rel4) <- c(N_age, N_strain)
  dim(n_Rvacc)  <- c(N_age, N_strain)
  dim(n_Rvacc_c) <- c(N_age, N_strain)

  


  
  
  ########## Transitions R1   ###########################################
  
  n_R1_wane[, ] <- Binomial(R1[i, j] - n_ageoR1[i, j] - n_muR1[i, j], p_RS2)
  
  n_R1vacc_c[, ] <- Binomial(R1[i, j] - n_ageoR1[i, j] - n_muR1[i, j]-  n_R1_wane[i,j],p_vacc[i]* switch_vacc3[i])
  n_R1vacc[, ]   <- n_R1vacc_c[i,j] * vacc_trans
  

    n_R1_tot[, ,] <- Binomial(((R1[i, j] - n_ageoR1[i, j] - n_muR1[i, j] - n_R1_wane[i, j]-n_R1vacc[i,j ]) ), 
  p_SE[i,k]*cross_map[j,k])
  
  
  # R1 reinfections from current strain
  #n_reinf_A2[,1]<- n_R1_rel1[i,1]
  n_R1_A1[,1]   <- n_R1_rel1[i,2]
  n_R1_A1[,2]   <- n_R1_rel2[i,3]
  n_R1_A1[,3]   <- n_R1_rel3[i,4]
  

  # R[1]
  
  n_R1_rel1[, 1] <- Binomial(sum(n_R1_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R1_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R1_tot[i, 1,]) - sum(n_R1_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R1_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R1_tot[i, 1,]) - sum(n_R1_rel1[i, 1:3])
  
  n_R1_E12[,1]  <- n_R1_rel1[i,3]
  n_R1_E12[,2]  <- n_R1_rel1[i,4]
  
  
  # R1[2]
  n_R1_rel2[, 1] <- Binomial(sum(n_R1_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R1_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R1_tot[i, 2,]) - sum(n_R1_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R1_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R1_tot[i, 2,]) - sum(n_R1_rel2[i, 1:3])
  
  
  #n_reinf_A3[,1]<- n_R1_rel2[i,1]
  n_R1_E13[, 1] <- n_R1_rel2[i, 2]
  n_R1_E13[, 2] <- n_R1_rel2[i, 4]
  
  
  # R1[3]
  n_R1_rel3[, 1] <- Binomial(sum(n_R1_tot[i, 3,]), rel_foi_strain[i, 1])
  n_R1_rel3[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R1_tot[i, 3,]) - sum(n_R1_rel3[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R1_rel3[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R1_tot[i, 3,]) - sum(n_R1_rel3[i, 1:3])
  
  #n_reinf_A4[,1]<- n_R1_rel3[i,1]
  n_R1_E14[, 1] <- n_R1_rel3[i, 2]
  n_R1_E14[, 2] <- n_R1_rel3[i, 3]
  
  
  
  dim(n_R1_tot)  <- c(N_age, N_strain_1,N_strain)
  dim(n_R1_rel1) <- c(N_age, N_strain)
  dim(n_R1_rel2) <- c(N_age, N_strain)
  dim(n_R1_rel3) <- c(N_age, N_strain)
  dim(n_R1vacc)  <- c(N_age, N_strain_1)
  dim(n_R1vacc_c) <- c(N_age, N_strain_1)
  ########## Transitions R2 
  
  n_R2_wane[, ] <- Binomial(R2[i, j] - n_ageoR2[i, j] - n_muR2[i, j], p_RS2)
  
  n_R2vacc_c[,] <- Binomial(R2[i, j] - n_ageoR2[i, j] - n_muR2[i, j]-n_R2_wane[i,j],p_vacc[i]* switch_vacc3[i])
  n_R2vacc[,] <- n_R2vacc_c[i,j] * vacc_trans
  
   n_R2_tot[, ,] <- Binomial(((R2[i, j] - n_ageoR2[i, j] - n_muR2[i, j] - n_R2_wane[i, j]-  n_R2vacc[i,j]) ), 
  p_SE[i,k]*cross_map[j,k])
  
  # R2 reinfections from current strain
  n_R2_A2[,1]   <- n_R2_rel1[i,1]
  n_R2_A2[,2]   <- n_R2_rel2[i,3]
  n_R2_A2[,3]   <- n_R2_rel3[i,4] 
  
  
  ## R2[1]
  n_R2_rel1[, 1] <- Binomial(sum(n_R2_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R2_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R2_tot[i, 1,]) - sum(n_R2_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R2_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R2_tot[i, 1,]) - sum(n_R2_rel1[i, 1:3])
  
  
  #n_reinf_A1[,1]<- n_R2_rel1[i,2]
  n_R2_E12[,1]  <- n_R2_rel1[i,3]
  n_R2_E12[,2]  <- n_R2_rel1[i,4]

  
  
  
  ## R2[2]
  n_R2_rel2[, 1] <- Binomial(sum(n_R2_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R2_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R2_tot[i, 2,]) - sum(n_R2_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R2_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R2_tot[i, 2,]) - sum(n_R2_rel2[i, 1:3])
  
  n_R2_E23[, 1] <- n_R2_rel2[i, 1]
  n_R2_E23[, 2] <- n_R2_rel2[i, 4]
  
  #n_reinf_A3[,2]<- n_R2_rel2[i,2]
  
  
  ## R2[3]
  n_R2_rel3[, 1] <- Binomial(sum(n_R2_tot[i, 3,]), rel_foi_strain[i, 1])
  n_R2_rel3[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R2_tot[i, 3,]) - sum(n_R2_rel3[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R2_rel3[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R2_tot[i, 3,]) - sum(n_R2_rel3[i, 1:3])
  
  
  n_R2_E24[, 1] <- n_R2_rel3[i, 1]
  n_R2_E24[, 2] <- n_R2_rel3[i, 3]

  
  #n_reinf_A4[,2]<-n_R2_rel3[i,2]
  
  dim(n_R2_tot)  <- c(N_age, N_strain_1, N_strain)
  dim(n_R2_rel1) <- c(N_age, N_strain)
  dim(n_R2_rel2) <- c(N_age, N_strain)
  dim(n_R2_rel3) <- c(N_age, N_strain)
  dim(n_R2vacc)  <- c(N_age,  N_strain_1)
  dim(n_R2vacc_c)  <- c(N_age,  N_strain_1)
  
  
  
  ########## Transitions R3 
  
  n_R3_wane[, ] <- Binomial(R3[i, j] - n_ageoR3[i, j] - n_muR3[i, j], p_RS2)
  
  n_R3vacc_c[, ]  <- Binomial(R3[i, j] - n_ageoR3[i, j] - n_muR3[i, j]-n_R3_wane[i,j],p_vacc[i]* switch_vacc3[i])
  n_R3vacc[, ]    <- n_R3vacc_c[i,j] * vacc_trans
  
  
  n_R3_tot[, , ] <- Binomial(((R3[i, j] - n_ageoR3[i, j] - n_muR3[i, j] - n_R3_wane[i, j]-n_R3vacc[i,j])), 
  p_SE[i,k]*cross_map[j,k])
  
  # R3 reinfections from current strain
  n_R3_A3[,1]   <- n_R3_rel1[i,1]
  n_R3_A3[,2]   <- n_R3_rel2[i,2]
  n_R3_A3[,3]   <- n_R3_rel3[i,4]
  
  # R3[1]
  n_R3_rel1[, 1] <- Binomial(sum(n_R3_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R3_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R3_tot[i, 1,]) - sum(n_R3_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R3_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R3_tot[i, 1,]) - sum(n_R3_rel1[i, 1:3])

   
  n_R3_E13[, 1] <- n_R3_rel1[i, 2]
  n_R3_E13[, 2] <- n_R3_rel1[i, 4]
  
  
  # R3[2]
  n_R3_rel2[, 1] <- Binomial(sum(n_R3_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R3_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R3_tot[i, 2,]) - sum(n_R3_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
   n_R3_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R3_tot[i, 2,]) - sum(n_R3_rel2[i, 1:3])
  
  n_R3_E23[, 1] <- n_R3_rel2[i, 1]
  n_R3_E23[, 2] <- n_R3_rel2[i, 4]
  
  
  # R3[3]
  n_R3_rel3[, 1] <- Binomial(sum(n_R3_tot[i, 3,]), rel_foi_strain[i, 1])
  n_R3_rel3[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R3_tot[i, 3,]) - sum(n_R3_rel3[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R3_rel3[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R3_tot[i, 3,]) - sum(n_R3_rel3[i, 1:3])

  
  n_R3_E34[,1]  <- n_R3_rel3[i,1]
  n_R3_E34[,2]  <- n_R3_rel3[i,2]
  
  
  dim(n_R3_tot) <- c(N_age, N_strain_1, N_strain)
  dim(n_R3_rel1) <- c(N_age, N_strain)
  dim(n_R3_rel2) <- c(N_age, N_strain)
  dim(n_R3_rel3) <- c(N_age, N_strain)
  dim(n_R3vacc)  <- c(N_age, N_strain_1)
  dim(n_R3vacc_c)  <- c(N_age, N_strain_1)
  ########## Transitions R4 
  
  n_R4_wane[, ] <- Binomial(R4[i, j] - n_ageoR4[i, j] - n_muR4[i, j], p_RS2)
  
  n_R4vacc_c[, ] <- Binomial(R4[i, j] - n_ageoR4[i, j] - n_muR4[i, j]-n_R4_wane[i,j],p_vacc[i]* switch_vacc3[i])
  n_R4vacc[, ]   <- n_R4vacc_c[i,j] * vacc_trans
  
 n_R4_tot[, ,] <- Binomial(((R4[i, j] - n_ageoR4[i, j] - n_muR4[i, j] - n_R4_wane[i, j]-n_R4vacc[i,j] ) ), 
  p_SE[i,k]*cross_map[j,k])
  
  
  # R4 reinfections from current strain
  n_R4_A4[,1]   <- n_R4_rel1[i,1]
  n_R4_A4[,2]   <- n_R4_rel2[i,2]
  n_R4_A4[,3]   <- n_R4_rel3[i,3]
  

  ## R4[1]
  n_R4_rel1[, 1] <- Binomial(sum(n_R4_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R4_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R4_tot[i, 1,]) - sum(n_R4_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R4_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R4_tot[i, 1,]) - sum(n_R4_rel1[i, 1:3])
 
  
  n_R4_E14[, 1] <- n_R4_rel1[i, 2]
  n_R4_E14[, 2] <- n_R4_rel1[i, 3]
  
  
  ## R4[2]
  n_R4_rel2[, 1] <- Binomial(sum(n_R4_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R4_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R4_tot[i, 2,]) - sum(n_R4_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R4_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R4_tot[i, 2,]) - sum(n_R4_rel2[i, 1:3])
  
  n_R4_E24[, 1] <- n_R4_rel2[i, 1]
  n_R4_E24[, 2] <- n_R4_rel2[i, 3]

  
  ## R4[3]
  n_R4_rel3[, 1] <- Binomial(sum(n_R4_tot[i, 3,]), rel_foi_strain[i, 1])
  n_R4_rel3[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R4_tot[i, 3,]) - sum(n_R4_rel3[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R4_rel3[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R4_tot[i, 3,]) - sum(n_R4_rel3[i, 1:3])

  
  n_R4_E34[,1]  <- n_R4_rel3[i,1]
  n_R4_E34[,2]  <- n_R4_rel3[i,2]
  

  
  dim(n_R4_tot) <- c(N_age, N_strain_1, N_strain)
  dim(n_R4_rel1) <- c(N_age, N_strain)
  dim(n_R4_rel2) <- c(N_age, N_strain)
  dim(n_R4_rel3) <- c(N_age, N_strain)
  dim(n_R4vacc)  <- c(N_age, N_strain_1)
  dim(n_R4vacc_c)  <- c(N_age, N_strain_1)
  
  
  ## TODO
  
  ########## Transitions R12 
  
  n_R12_wane[, ] <- Binomial(R12[i, j] - n_ageoR12[i, j] - n_muR12[i, j], p_RS2)
  
  n_R12vacc_c[, ] <- Binomial(R12[i, j] - n_ageoR12[i, j] - n_muR12[i, j]- n_R12_wane[i,j],p_vacc[i]* switch_vacc3[i])
  n_R12vacc[, ]   <- n_R12vacc_c[i,j] * vacc_trans
  
  n_R12_tot[, ,] <- Binomial(((R12[i, j] - n_ageoR12[i, j] - n_muR12[i, j] - n_R12_wane[i, j]- n_R12vacc[i, j]) ), 
  p_SE[i,k]*cross_map[j,k])
  
  # Refinfection with same strain
  n_R12_A12[,1]   <- n_R12_rel1[i,3]
  n_R12_A12[,2]   <- n_R12_rel2[i,4] 
  
  # R12[1]
  n_R12_rel1[, 1] <- Binomial(sum(n_R12_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R12_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R12_tot[i, 1,]) - sum(n_R12_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R12_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R12_tot[i, 1,]) - sum(n_R12_rel1[i, 1:3])

  
  
  
  
  # R12[2]
  n_R12_rel2[, 1] <- Binomial(sum(n_R12_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R12_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R12_tot[i, 2,]) - sum(n_R12_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R12_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R12_tot[i, 2,]) - sum(n_R12_rel2[i, 1:3])
  
  n_R12_E4rd[, 1] <-0
  n_R12_E4rd[, 2] <-0
  n_R12_E4rd[, 4] <- n_R12_rel1[i, 4]
  n_R12_E4rd[, 3] <- n_R12_rel2[i, 3]

  
  
  
  dim(n_R12_rel1) <- c(N_age, N_strain)
  dim(n_R12_rel2) <- c(N_age, N_strain)
  dim(n_R13_rel1) <- c(N_age, N_strain)
  dim(n_R13_rel2) <- c(N_age, N_strain)
  dim(n_R14_rel1) <- c(N_age, N_strain)
  dim(n_R14_rel2) <- c(N_age, N_strain)
  dim(n_R23_rel1) <- c(N_age, N_strain)
  dim(n_R23_rel2) <- c(N_age, N_strain)
  dim(n_R24_rel1) <- c(N_age, N_strain)
  dim(n_R24_rel2) <- c(N_age, N_strain)
  dim(n_R34_rel1) <- c(N_age, N_strain)
  dim(n_R34_rel2) <- c(N_age, N_strain)
  
  
  dim(n_R12_tot) <- c(N_age, N_strain_2, N_strain)
  dim(n_R13_tot) <- c(N_age, N_strain_2, N_strain)
  dim(n_R14_tot) <- c(N_age, N_strain_2, N_strain)
  dim(n_R23_tot) <- c(N_age, N_strain_2, N_strain)
  dim(n_R24_tot) <- c(N_age, N_strain_2, N_strain)
  dim(n_R34_tot) <- c(N_age, N_strain_2, N_strain)
  
  dim(n_R12vacc) <- c(N_age, N_strain_2)
  dim(n_R13vacc) <- c(N_age, N_strain_2)
  dim(n_R14vacc) <- c(N_age, N_strain_2)
  dim(n_R23vacc) <- c(N_age, N_strain_2)
  dim(n_R24vacc) <- c(N_age, N_strain_2)
  dim(n_R34vacc)<-  c(N_age,N_strain_2)
  
  dim(n_R12vacc_c) <- c(N_age, N_strain_2)
  dim(n_R13vacc_c) <- c(N_age, N_strain_2)
  dim(n_R14vacc_c) <- c(N_age, N_strain_2)
  dim(n_R23vacc_c) <- c(N_age, N_strain_2)
  dim(n_R24vacc_c) <- c(N_age, N_strain_2)
  dim(n_R34vacc_c)<-  c(N_age,N_strain_2)
  
  
  
  ########## Transitions R13 
  
  n_R13_wane[, ] <- Binomial(R13[i, j] - n_ageoR13[i, j] - n_muR13[i, j], p_RS2)
  
  n_R13vacc_c[, ] <- Binomial(R13[i, j] - n_ageoR13[i, j] - n_muR13[i, j]-  n_R13_wane[i,j], p_RS2) 
  
  n_R13vacc[, ] <- n_R13vacc_c[i,j] * vacc_trans
  
   n_R13_tot[, ,] <- Binomial(((R13[i, j] - n_ageoR13[i, j] - n_muR13[i, j] - n_R13_wane[i, j]-n_R13vacc[i,j ] ) ), 
  p_SE[i,k]*cross_map[j,k])
  
    # reinfection from same strain
  n_R13_A13[,1]   <- n_R13_rel1[i,2]
  n_R13_A13[,2]   <- n_R13_rel2[i,4]
  # R13[1]
  n_R13_rel1[, 1] <- Binomial(sum(n_R13_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R13_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R13_tot[i, 1,]) - sum(n_R13_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R13_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R13_tot[i, 1,]) - sum(n_R13_rel1[i, 1:3])
 
  

  
  
  
  # R13[2]
  n_R13_rel2[, 1] <- Binomial(sum(n_R13_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R13_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R13_tot[i, 2,]) - sum(n_R13_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R13_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R13_tot[i, 2,]) - sum(n_R13_rel2[i, 1:3])
 
  n_R13_E4rd[, 1] <- 0
  n_R13_E4rd[, 2] <- n_R13_rel2[i, 2]
  n_R13_E4rd[, 3] <- 0
  n_R13_E4rd[, 4] <- n_R13_rel1[i, 4]  
  
  ########## Transitions R14 (not 1 and 4, j 2 and 3) 
  
  n_R14_wane[, ] <- Binomial(R14[i, j] - n_ageoR14[i, j] - n_muR14[i, j], p_RS2)
  
  n_R14vacc_c[, ] <- Binomial(R14[i, j] - n_ageoR14[i, j] - n_muR14[i, j]-n_R14_wane[i,j], p_RS2)
  n_R14vacc[, ]   <- n_R14vacc_c[i,j] * vacc_trans
  
  
  n_R14_tot[, ,] <- Binomial(((R14[i, j] - n_ageoR14[i, j] - n_muR14[i, j] - n_R14_wane[i, j]- n_R14vacc[i,j ]  ) ), 
  p_SE[i,k]*cross_map[j,k]) 

   # reinfection from same strain
  n_R14_A14[,1]   <- n_R14_rel1[i,2]
  n_R14_A14[,2]   <- n_R14_rel2[i,3] 
  
  
  # R14[1]
  n_R14_rel1[, 1] <- Binomial(sum(n_R14_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R14_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R14_tot[i, 1,]) - sum(n_R14_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R14_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R14_tot[i, 1,]) - sum(n_R14_rel1[i, 1:3])

  
  
  
  
  # R14[2]
  n_R14_rel2[, 1] <- Binomial(sum(n_R14_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R14_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R14_tot[i, 2,]) - sum(n_R14_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R14_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R14_tot[i, 2,]) - sum(n_R14_rel2[i, 1:3])

  n_R14_E4rd[, 1] <- 0
  n_R14_E4rd[, 2] <- n_R14_rel2[i, 2]
  n_R14_E4rd[, 3] <- n_R14_rel1[i, 3]
  n_R14_E4rd[, 4] <- 0
  
  
  ########## Transitions R23 (not 2 and 3, j 1 and 4) 
  
  n_R23_wane[, ] <- Binomial(R23[i, j] - n_ageoR23[i, j] - n_muR23[i, j], p_RS2)
  
  n_R23vacc_c[, ] <- Binomial(R23[i, j] - n_ageoR23[i, j] - n_muR23[i, j]-  n_R23_wane[i,j], p_RS2)
  n_R23vacc[, ]   <- n_R23vacc_c[i,j]*vacc_trans
  
  
  n_R23_tot[, ,] <- Binomial(((R23[i, j] - n_ageoR23[i, j] - n_muR23[i, j] - n_R23_wane[i, j]) ), 
  p_SE[i,k]*cross_map[j,k])

  # reinf from same strain
  n_R23_A23[,1]   <- n_R23_rel1[i,1]
  n_R23_A23[,2]   <- n_R23_rel2[i,4]
  
  # R23[1]
  n_R23_rel1[, 1] <- Binomial(sum(n_R23_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R23_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R23_tot[i, 1,]) - sum(n_R23_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R23_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R23_tot[i, 1,]) - sum(n_R23_rel1[i, 1:3])
   
  
  
  
  # R23[2]
  n_R23_rel2[, 1] <- Binomial(sum(n_R23_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R23_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R23_tot[i, 2,]) - sum(n_R23_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R23_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R23_tot[i, 2,]) - sum(n_R23_rel2[i, 1:3])
  
  
   n_R23_E4rd[, 1] <- n_R23_rel2[i, 1]
   n_R23_E4rd[, 2] <- 0
   n_R23_E4rd[, 3] <- 0
   n_R23_E4rd[, 4] <- n_R23_rel1[i, 4]
  
  ########## Transitions R24 (not 2 and 4, j 1 and 3) 
  
  n_R24_wane[, ] <- Binomial(R24[i, j] - n_ageoR24[i, j] - n_muR24[i, j], p_RS2)
  
  n_R24vacc_c[, ] <- Binomial(R24[i, j] - n_ageoR24[i, j] - n_muR24[i, j]-n_R24_wane[i,j], p_RS2)
  n_R24vacc[, ]   <- n_R24vacc_c[i,j]*vacc_trans
  
  n_R24_tot[, ,] <- Binomial(((R24[i, j] - n_ageoR24[i, j] - n_muR24[i, j] - n_R24_wane[i, j]-n_R24vacc[i, j] ) ), 
  p_SE[i,k]*cross_map[j,k])
  
  # Reinfection with same strain
  n_R24_A24[,1]   <- n_R24_rel1[i,1]
  n_R24_A24[,2]   <- n_R24_rel2[i,3]
  
  
  # R24[1]
  n_R24_rel1[, 1] <- Binomial(sum(n_R24_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R24_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R24_tot[i, 1,]) - sum(n_R24_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R24_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R24_tot[i, 1,]) - sum(n_R24_rel1[i, 1:3]) 


  
  
  
  # R24[2]
  n_R24_rel2[, 1] <- Binomial(sum(n_R24_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R24_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R24_tot[i, 2,]) - sum(n_R24_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R24_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R24_tot[i, 2,]) - sum(n_R24_rel2[i, 1:3])
  
  n_R24_E4rd[, 1] <- n_R24_rel2[i, 1]
  n_R24_E4rd[, 2] <- 0
  n_R24_E4rd[, 3] <- n_R24_rel1[i, 3]  
  n_R24_E4rd[, 4] <- 0
  
  ########## Transitions R34 (not 3 and 4, j 1 and 2) 
  
  n_R34_wane[, ] <- Binomial(R34[i, j] - n_ageoR34[i, j] - n_muR34[i, j], p_RS2)
  
  n_R34vacc_c[, ] <- Binomial(R34[i, j] - n_ageoR34[i, j] - n_muR34[i, j]-n_R34_wane[i,j], p_RS2)
  n_R34vacc[, ]   <- n_R34vacc_c[i,j]*vacc_trans
  
  
  n_R34_tot[, ,] <- Binomial(((R34[i, j] - n_ageoR34[i, j] - n_muR34[i, j] - n_R34_wane[i, j]-n_R34vacc[i, j] ) ), 
  p_SE[i,k]*cross_map[j,k])
  
  # reinf from strain
  n_R34_A34[,1]   <- n_R34_rel1[i,1]
  n_R34_A34[,2]   <- n_R34_rel2[i,2]
  
  # R34[1]
  n_R34_rel1[, 1] <- Binomial(sum(n_R34_tot[i, 1,]), rel_foi_strain[i, 1])
  n_R34_rel1[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R34_tot[i, 1,]) - sum(n_R34_rel1[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R34_rel1[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R34_tot[i, 1,]) - sum(n_R34_rel1[i, 1:3])   

  
  
  
  # R34[2]
  n_R34_rel2[, 1] <- Binomial(sum(n_R34_tot[i, 2,]), rel_foi_strain[i, 1])
  n_R34_rel2[, 2:3] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    Binomial(sum(n_R34_tot[i, 2,]) - sum(n_R34_rel2[i, 1:(j - 1)]), rel_foi_strain[i, j]/sum(rel_foi_strain[i, j:N_strain]))
  n_R34_rel2[, 4] <- if (sum(rel_foi_strain[i, j:N_strain]) == 0) 0 else 
    sum(n_R34_tot[i, 2,]) - sum(n_R34_rel2[i, 1:3])
  
  n_R34_E4rd[, 1] <- n_R34_rel2[i, 1] 
  n_R34_E4rd[, 2] <- n_R34_rel1[i, 2]
  n_R34_E4rd[, 3] <- 0
  n_R34_E4rd[, 4] <- 0
  
  ########### R 4rd transitions HERE
  
  n_R4rd_wane[, ] <- Binomial(R4rd[i, j] - n_ageoR4rd[i, j] - n_muR4rd[i, j], p_RS2)
  
  n_R4rdvacc_c[, ] <- Binomial(R4rd[i, j] - n_ageoR4rd[i, j] - n_muR4rd[i, j]-  n_R4rd_wane[i,j], p_RS2)  
  n_R4rdvacc[, ]   <- n_R4rdvacc_c[i,j] * vacc_trans  
  
  # Leave R 
  n_R4rd_tot[, ,] <- Binomial(((R4rd[i, j] -
   n_ageoR4rd[i, j] - 
   n_muR4rd[i, j] - 
   n_R4rd_wane[i, j]-
   n_R4rdvacc[i, j] ) ), p_SE[i,k]*cross_map[j,k])
  
  # reinfections from same strain
  n_R4rd_A4rd[, ] <-  sum(n_R4rd_tot[i,,j])  


  dim(n_R4rd_tot) <- c(N_age, N_strain, N_strain)
  dim(n_R4rdvacc) <- c(N_age, N_strain)
  dim(n_R4rdvacc_c) <- c(N_age, N_strain)
  
  
  ############ Reinfection arrays together
  n_reinf_A1[,1]<- n_R2_rel1[i,2]
  n_reinf_A1[,2]<- n_R3_rel1[i,3]
  n_reinf_A1[,3]<- n_R4_rel1[i,4]
  
  n_reinf_A2[,1]<- n_R1_rel1[i,1]
  n_reinf_A2[,2]<- n_R3_rel2[i,3]
  n_reinf_A2[,3]<- n_R4_rel2[i,4]
  
  n_reinf_A3[,1]<- n_R1_rel2[i,1]
  n_reinf_A3[,2]<- n_R2_rel2[i,2]
  n_reinf_A3[,3]<- n_R4_rel3[i,4]
  
  n_reinf_A4[,1]<- n_R1_rel3[i,1]
  n_reinf_A4[,2]<-n_R2_rel3[i,2]
  n_reinf_A4[,3]<- n_R3_rel3[i,3]
  
  n_reinf_A12[,1] <- n_R23_rel1[i,3]+n_R13_rel1[i,3]
  n_reinf_A12[,2] <- n_R24_rel1[i,4]+n_R14_rel1[i,4]
  
  
  n_reinf_A13[,1] <- n_R23_rel1[i,2]+n_R12_rel1[i,2]
  n_reinf_A13[,2] <- n_R14_rel2[i,4]+n_R34_rel1[i,4]
  
  n_reinf_A14[,1] <- n_R12_rel2[i,2]+n_R24_rel1[i,2]
  n_reinf_A14[,2] <- n_R34_rel1[i,3]+n_R13_rel2[i,3]
  
  n_reinf_A23[,1] <- n_R13_rel1[i,1]+n_R12_rel1[i,1]
  n_reinf_A23[,2] <- n_R34_rel2[i,4]+n_R24_rel2[i,4]
  
  n_reinf_A24[,1] <- n_R14_rel1[i,1]+n_R12_rel2[i,1]
  n_reinf_A24[,2] <- n_R34_rel2[i,3]+n_R23_rel2[i,3]
  
  
  n_reinf_A34[,1] <- n_R14_rel2[i,1]+n_R13_rel2[i,1]
  n_reinf_A34[,2] <- n_R24_rel2[i,2]+n_R23_rel2[i,2]
  
  
  
 
  
  
  ##### Deaths 
  ###################################
  
  n_allDeath<- 
    sum(n_muG)+
    sum(n_muM)+
    sum(n_muS)+
    sum(n_muE)+
    sum(n_muI)+
    sum(n_muA)+
    sum(n_muR)+
    sum(n_muE1)+
    sum(n_muI1)+
    sum(n_muA1)+
    sum(n_muR1)+
    sum(n_muE2)+
    sum(n_muI2)+
    sum(n_muA2)+
    sum(n_muR2)+
    sum(n_muE3)+
    sum(n_muI3)+
    sum(n_muA3)+
    sum(n_muR3)+
    sum(n_muE4)+
    sum(n_muI4)+
    sum(n_muA4)+
    sum(n_muR4)+
    sum(n_muE12)+
    sum(n_muI12)+
    sum(n_muA12)+
    sum(n_muR12)+
    sum(n_muE13)+
    sum(n_muI13)+
    sum(n_muA13)+
    sum(n_muR13)+
    sum(n_muE14)+
    sum(n_muI14)+
    sum(n_muA14)+
    sum(n_muR14)+
    sum(n_muE23)+
    sum(n_muI23)+
    sum(n_muA23)+
    sum(n_muR23)+
    sum(n_muE24)+
    sum(n_muI24)+
    sum(n_muA24)+
    sum(n_muR24)+
    sum(n_muE34)+
    sum(n_muI34)+
    sum(n_muA34)+
    sum(n_muR34)+
    sum(n_muE4rd)+
    sum(n_muI4rd)+
    sum(n_muA4rd)+
    sum(n_muR4rd)
  
  
  # Births to keep stable population equal to deaths
  n_bG[] <- if (i==1) (n_allDeath*p_nonsecretor) else 0 
  n_bM[] <- if (i==1) (n_allDeath*(1-p_nonsecretor)*prev) else 0 
  n_bS[] <- if (i==1) ((n_allDeath-n_bG[1]-n_bM[1])) else 0 
  
  
  
  # Vaccination transitions - see index_states excel  -----------------------
  # Transitions imply a vaccine for GI3 and GIi4 and states are allocated accordignly
  
  
  n_vacR1[,1]   <- 0
  n_vacR1[,2]   <- n_Mvacc[i] + n_Svacc[i] + n_Rvacc[i,1] + n_R1vacc[i,2]
  n_vacR1[,3]   <- 0
  
  n_vacR3[,1]   <- n_Rvacc[i,3] + n_R3vacc[i,1]
  n_vacR3[,2]   <- 0 
  
  
  n_vacR12[,1]  <- n_R1vacc[i,1] + n_R2vacc[i,1] + n_R12vacc[i,1]
  n_vacR12[,2]  <- 0
  
  n_vacR13[,1]  <- n_Rvacc[i,2] + n_R13vacc[i,1] 
  n_vacR13[,2]  <- n_Rvacc[i,4] + n_R13vacc[i,2] 
  
  n_vacR14[,1]  <- 0
  n_vacR14[,2]  <- n_R1vacc[i,3] + n_R4vacc[i,1] + n_R14vacc[i,2]
  
  n_vacR23[,1]  <- n_R2vacc[i,2] + n_R3vacc[i,2] + n_R23vacc[i,1]
  n_vacR23[,2]  <- 0
  
  n_vacR34[,1]   <- n_R3vacc[i,3] + n_R4vacc[i,3] + n_R34vacc[i,1]
  n_vacR34[,2]   <- 0 
  
  
  n_vacR1234[,1]<- n_R2vacc[i,3]+n_R4vacc[i,2]+n_R12vacc[i,2]+n_R23vacc[i,2]+n_R24vacc[i,2]+n_R34vacc[i,2]
  n_vacR1234[,2]<- 0
  n_vacR1234[,3]<- n_R14vacc[i,1]+n_R24vacc[i,1]
  n_vacR1234[,4]<- 0
  
  
  # E to A
  n_vacA1[,1]   <- 0
  n_vacA1[,2]   <- n_Evacc[i,1] + n_E1vacc[i,2]
  n_vacA1[,3]   <- 0
  
  n_vacA3[,1]   <- n_Evacc[i,3] + n_E3vacc[i,1]
  n_vacA3[,2]   <- 0 
  
  n_vacA12[,1]  <- n_E1vacc[i,1] + n_E2vacc[i,1] + n_E12vacc[i,1]
  n_vacA12[,2]  <- 0
  
  n_vacA13[,1]  <- n_Evacc[i,2] + n_E13vacc[i,1] 
  n_vacA13[,2]  <- n_Evacc[i,4] + n_E13vacc[i,2] 
  
  n_vacA14[,1]  <- 0
  n_vacA14[,2]  <- n_E1vacc[i,3] + n_E4vacc[i,1] + n_E14vacc[i,2]
  
  n_vacA23[,1]  <- n_E2vacc[i,2] + n_E3vacc[i,2] + n_E23vacc[i,1]
  n_vacA23[,2]  <- 0
  
  n_vacA34[,1]   <- n_E3vacc[i,3] + n_E4vacc[i,3] + n_E34vacc[i,1]
  n_vacA34[,2]   <- 0 
  
  n_vacA1234[,1]<- n_E2vacc[i,3]+n_E4vacc[i,2]+n_E12vacc[i,2]+n_E23vacc[i,2]+n_E24vacc[i,2]+n_E34vacc[i,2]
  n_vacA1234[,2]<- 0
  n_vacA1234[,3]<- n_E14vacc[i,1]+n_E24vacc[i,1]
  n_vacA1234[,4]<- 0
  
  
  
  #Doses
  update(n_vacc_doses)<- 
    sum(n_Svacc_c)+
    sum(n_Mvacc_c)+
    sum(n_Rvacc_c)+
    sum(n_R1vacc_c)+
    sum(n_R2vacc_c)+
    sum(n_R3vacc_c)+
    sum(n_R4vacc_c)+
    sum(n_R12vacc_c)+
    sum(n_R13vacc_c)+
    sum(n_R14vacc_c)+
    sum(n_R23vacc_c)+
    sum(n_R24vacc_c)+
    sum(n_R34vacc_c)+
    sum(n_R4rdvacc_c)+
    sum(n_Evacc	)+
    sum(n_E1vacc)+
    sum(n_E2vacc)+
    sum(n_E3vacc)+
    sum(n_E4vacc)+
    sum(n_E12vacc)+
    sum(n_E13vacc)+
    sum(n_E14vacc)+
    sum(n_E23vacc)+
    sum(n_E24vacc)+
    sum(n_E34vacc)+
    sum(n_E4rdvacc)
  
  
  
  
  dim(n_vacR1)<-c(N_age,N_strain_1)
  dim(n_vacR3)<-c(N_age,N_strain_1)
  dim(n_vacR12)<-c(N_age,N_strain_2)
  dim(n_vacR1234)<-c(N_age,N_strain)
  dim(n_vacR13)<-c(N_age,N_strain_2)
  dim(n_vacR14)<-c(N_age,N_strain_2)
  dim(n_vacR23)<-c(N_age,N_strain_2)
  dim(n_vacR34)<-c(N_age,N_strain_2)
  
  dim(n_vacA1)<-c(N_age,N_strain_1)
  dim(n_vacA3)<-c(N_age,N_strain_1)
  dim(n_vacA12)<-c(N_age,N_strain_2)
  dim(n_vacA1234)<-c(N_age,N_strain)
  dim(n_vacA13)<-c(N_age,N_strain_2)
  dim(n_vacA14)<-c(N_age,N_strain_2)
  dim(n_vacA23)<-c(N_age,N_strain_2)
  dim(n_vacA34)<-c(N_age,N_strain_2)
  
  
  
  ## Initial states:
  initial(G[])<-pop[i]*p_nonsecretor
  initial(M[])<-0
  initial(S[])<-pop[i]*(1-p_nonsecretor)
  initial(E[,])<-0
  initial(I[,])<-if(i==5) 10 else 0
  initial(A[,])<-0
  initial(R[,])<-0
  initial(E1[,])<-0
  initial(I1[,])<-0
  initial(A1[,])<-0
  initial(R1[,])<-0
  initial(E2[,])<-0
  initial(I2[,])<-0
  initial(A2[,])<-0
  initial(R2[,])<-0
  initial(E3[,])<-0
  initial(I3[,])<-0
  initial(A3[,])<-0
  initial(R3[,])<-0
  initial(E4[,])<-0
  initial(I4[,])<-0
  initial(A4[,])<-0
  initial(R4[,])<-0
  initial(E12[,])<-0
  initial(I12[,])<-0
  initial(A12[,])<-0
  initial(R12[,])<-0
  initial(E13[,])<-0
  initial(I13[,])<-0
  initial(A13[,])<-0
  initial(R13[,])<-0
  initial(E14[,])<-0
  initial(I14[,])<-0
  initial(A14[,])<-0
  initial(R14[,])<-0
  initial(E23[,])<-0
  initial(I23[,])<-0
  initial(A23[,])<-0
  initial(R23[,])<-0
  initial(E24[,])<-0
  initial(I24[,])<-0
  initial(A24[,])<-0
  initial(R24[,])<-0
  initial(E34[,])<-0
  initial(I34[,])<-0
  initial(A34[,])<-0
  initial(R34[,])<-0
  initial(E4rd[,])<-0
  initial(I4rd[,])<-0
  initial(A4rd[,])<-0
  initial(R4rd[,])<-0
  
  
  # Outputs
  initial(inc_day_gii4[]) <- 0
  initial(inc_day_gii[]) <- 0
  initial(inc_day_gi3[]) <- 0
  initial(inc_day_gi[]) <- 0
  initial(inc_day_all_gii4[]) <- 0
  initial(inc_day_all_gii[]) <- 0
  initial(inc_day_all_gi3[]) <- 0
  initial(inc_day_all_gi[]) <- 0
  initial(inc_day_all[]) <- 0
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
  initial(death_reported_wk, zero_every = 7) <- 0
  initial(death_reported_wk_all[], zero_every = 7) <- 0
  initial(reported_wk_gi3, zero_every = 7) <- 0
  initial(reported_wk_gi, zero_every = 7) <- 0
  initial(reported_wk_gii4, zero_every = 7) <- 0
  initial(reported_wk_gii, zero_every = 7) <- 0
  initial(covid_timeline)<-0
  initial(new_cases)<-0
  initial(new_hosp_adult_yr, zero_every = 365) <- 0
  initial(new_hosp_elder_yr, zero_every = 365) <- 0
  initial(new_cases_week, zero_every = 7) <- 0
  initial(new_cases_week_gi, zero_every = 7) <- 0
  initial(new_cases_week_gii, zero_every = 7) <- 0
  initial(new_cases_week_gi3, zero_every = 7) <- 0
  initial(new_cases_week_gii4, zero_every = 7) <- 0
  initial(seasonality)<-0
  initial(sus_gi)<-0
  initial(sus_gii)<-0
  initial(sus_gi3)<-0
  initial(sus_gii4)<-0
  initial(n_vacc_doses)<-0
  
  initial(sus_cross_gi)<-0
  initial(sus_cross_gii)<-0
  initial(sus_cross_gi3)<-0
  initial(sus_cross_gii4)<-0
  
  initial(sus_re_gi)<-0
  initial(sus_re_gii)<-0
  initial(sus_re_gi3)<-0
  initial(sus_re_gii4)<-0
  initial(pop_all[])<-0
  initial(foi_gi3[])<-0
  initial(foi_gi[])<-0
  initial(foi_gii4[])<-0
  initial(foi_gii[])<-0
  
  
  
  
  # Parameters --------------------------------------------------------------
  
  
  # Calibrated Parameters --------------------------------------------------------
  beta_1 <- parameter(0.2)
  beta_2 <- parameter(0.2)
  beta_3 <- parameter(0.2)
  beta_4 <- parameter(0.2)
  maternalAB <- parameter(180)
  aduRR <- parameter(0.1)
  imm_yr <- parameter(5.1)
  imm_fac <- parameter(2)
  repfac_0 <- parameter(287)
  repfac_5 <- parameter(287)
  repfac_15 <- parameter(287)
  repfac_65p <- parameter(287)
  crossp_GI <-  parameter(0.05)
  crossp_GII <- parameter(0.05)
  season_lag  <-parameter(30)
  season_amp  <-parameter(0.15)
  
  
  # Fixed with default--------------------------------------------------------
  epsilon <- parameter(1)
  theta_5 <- parameter(2.5)
  theta_5p <- parameter(1.5)
  sigma <- parameter(15)
  rr_inf_asymp <- parameter(0.05)
  p_nonsecretor <- parameter(0.2)
  #sus_vec[]       <- 1#if (i==3) 1.0-p_nonsecretor else as.numeric(1)
  #sus_vecnew[]       <- if (i==3) 1.0-p_nonsecretor else as.numeric(1)
  geno_frac       <-0.2
  N_strain        <- parameter(4)
  t_zero          <-parameter(0)
  t_campaign_stop <-parameter(180)
  hosp_rate_adult <-parameter(4.8/1e4)
  hosp_rate_elder <-parameter(43/1e4)
  
  # Fixed with NO default--------------------------------------------------------
  pop<- parameter()
  mu <- parameter()
  cfr<- parameter()
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
  vaccination_coverage<- parameter()
  vacc_switch_on      <- parameter()
  campaign_switch     <- parameter()
  vacc_trans          <- parameter()
  vacc_dis            <- parameter()  
  
  
  
  # Dimensions--------------------------------------------------------
  N_strain_1<-N_strain-1
  N_strain_2<-N_strain-2
  dim(aging_vec)<- N_age
    dim(	G	)<-N_age

  dim(vaccination_coverage)<- N_age
  dim(vacc_switch_on)      <- N_age
  dim(campaign_switch)     <- N_age  
  dim(	M	)<-N_age
  dim(	S	)<-N_age
  dim(	E	)<-c(N_age,N_strain)
  dim(	I	)<-c(N_age,N_strain)
  dim(	A	)<-c(N_age,N_strain)
  dim(	R	)<-c(N_age,N_strain)
  dim(	E1	)<-c(N_age,N_strain_1)
  dim(	I1	)<-c(N_age,N_strain_1)
  dim(	A1	)<-c(N_age,N_strain_1)
  dim(	R1	)<-c(N_age,N_strain_1)
  dim(	E2	)<-c(N_age,N_strain_1)
  dim(	I2	)<-c(N_age,N_strain_1)
  dim(	A2	)<-c(N_age,N_strain_1)
  dim(	R2	)<-c(N_age,N_strain_1)
  dim(	E3	)<-c(N_age,N_strain_1)
  dim(	I3	)<-c(N_age,N_strain_1)
  dim(	A3	)<-c(N_age,N_strain_1)
  dim(	R3	)<-c(N_age,N_strain_1)
  dim(	E4	)<-c(N_age,N_strain_1)
  dim(	I4	)<-c(N_age,N_strain_1)
  dim(	A4	)<-c(N_age,N_strain_1)
  dim(	R4	)<-c(N_age,N_strain_1)
  dim(	E12	)<-c(N_age,N_strain_2)
  dim(	I12	)<-c(N_age,N_strain_2)
  dim(	A12	)<-c(N_age,N_strain_2)
  dim(	R12	)<-c(N_age,N_strain_2)
  dim(	E13	)<-c(N_age,N_strain_2)
  dim(	I13	)<-c(N_age,N_strain_2)
  dim(	A13	)<-c(N_age,N_strain_2)
  dim(	R13	)<-c(N_age,N_strain_2)
  dim(	E14	)<-c(N_age,N_strain_2)
  dim(	I14	)<-c(N_age,N_strain_2)
  dim(	A14	)<-c(N_age,N_strain_2)
  dim(	R14	)<-c(N_age,N_strain_2)
  dim(	E23	)<-c(N_age,N_strain_2)
  dim(	I23	)<-c(N_age,N_strain_2)
  dim(	A23	)<-c(N_age,N_strain_2)
  dim(	R23	)<-c(N_age,N_strain_2)
  dim(	E24	)<-c(N_age,N_strain_2)
  dim(	I24	)<-c(N_age,N_strain_2)
  dim(	A24	)<-c(N_age,N_strain_2)
  dim(	R24	)<-c(N_age,N_strain_2)
  dim(	E34	)<-c(N_age,N_strain_2)
  dim(	I34	)<-c(N_age,N_strain_2)
  dim(	A34	)<-c(N_age,N_strain_2)
  dim(	R34	)<-c(N_age,N_strain_2)
  dim(	E4rd	)<-c(N_age,N_strain)
  dim(	I4rd	)<-c(N_age,N_strain)
  dim(	A4rd	)<-c(N_age,N_strain)
  dim(	R4rd	)<-c(N_age,N_strain)
  
  dim(n_ageiG)<-N_age
  dim(n_ageiM)<-N_age
  dim(n_ageiS)<-N_age
  dim(n_ageiE)<-c(N_age,N_strain)
  dim(n_ageiI)<-c(N_age,N_strain)
  dim(n_ageiA)<-c(N_age,N_strain)
  dim(n_ageiR)<-c(N_age,N_strain)
  dim(n_ageiE1)<-c(N_age,N_strain_1)
  dim(n_ageiI1)<-c(N_age,N_strain_1)
  dim(n_ageiA1)<-c(N_age,N_strain_1)
  dim(n_ageiR1)<-c(N_age,N_strain_1)
  dim(n_ageiE2)<-c(N_age,N_strain_1)
  dim(n_ageiI2)<-c(N_age,N_strain_1)
  dim(n_ageiA2)<-c(N_age,N_strain_1)
  dim(n_ageiR2)<-c(N_age,N_strain_1)
  dim(n_ageiE3)<-c(N_age,N_strain_1)
  dim(n_ageiI3)<-c(N_age,N_strain_1)
  dim(n_ageiA3)<-c(N_age,N_strain_1)
  dim(n_ageiR3)<-c(N_age,N_strain_1)
  dim(n_ageiE4)<-c(N_age,N_strain_1)
  dim(n_ageiI4)<-c(N_age,N_strain_1)
  dim(n_ageiA4)<-c(N_age,N_strain_1)
  dim(n_ageiR4)<-c(N_age,N_strain_1)
  dim(n_ageiE12)<-c(N_age,N_strain_2)
  dim(n_ageiI12)<-c(N_age,N_strain_2)
  dim(n_ageiA12)<-c(N_age,N_strain_2)
  dim(n_ageiR12)<-c(N_age,N_strain_2)
  dim(n_ageiE13)<-c(N_age,N_strain_2)
  dim(n_ageiI13)<-c(N_age,N_strain_2)
  dim(n_ageiA13)<-c(N_age,N_strain_2)
  dim(n_ageiR13)<-c(N_age,N_strain_2)
  dim(n_ageiE14)<-c(N_age,N_strain_2)
  dim(n_ageiI14)<-c(N_age,N_strain_2)
  dim(n_ageiA14)<-c(N_age,N_strain_2)
  dim(n_ageiR14)<-c(N_age,N_strain_2)
  dim(n_ageiE23)<-c(N_age,N_strain_2)
  dim(n_ageiI23)<-c(N_age,N_strain_2)
  dim(n_ageiA23)<-c(N_age,N_strain_2)
  dim(n_ageiR23)<-c(N_age,N_strain_2)
  dim(n_ageiE24)<-c(N_age,N_strain_2)
  dim(n_ageiI24)<-c(N_age,N_strain_2)
  dim(n_ageiA24)<-c(N_age,N_strain_2)
  dim(n_ageiR24)<-c(N_age,N_strain_2)
  dim(n_ageiE34)<-c(N_age,N_strain_2)
  dim(n_ageiI34)<-c(N_age,N_strain_2)
  dim(n_ageiA34)<-c(N_age,N_strain_2)
  dim(n_ageiR34)<-c(N_age,N_strain_2)
  dim(n_ageiE4rd)<-c(N_age,N_strain)
  dim(n_ageiI4rd)<-c(N_age,N_strain)
  dim(n_ageiA4rd)<-c(N_age,N_strain)
  dim(n_ageiR4rd)<-c(N_age,N_strain)
  
  dim(n_ageoG)<-N_age
  dim(n_ageoM)<-N_age
  dim(n_ageoS)<-N_age
  dim(n_ageoE)<-c(N_age,N_strain)
  dim(n_ageoI)<-c(N_age,N_strain)
  dim(n_ageoA)<-c(N_age,N_strain)
  dim(n_ageoR)<-c(N_age,N_strain)
  dim(n_ageoE1)<-c(N_age,N_strain_1)
  dim(n_ageoI1)<-c(N_age,N_strain_1)
  dim(n_ageoA1)<-c(N_age,N_strain_1)
  dim(n_ageoR1)<-c(N_age,N_strain_1)
  dim(n_ageoE2)<-c(N_age,N_strain_1)
  dim(n_ageoI2)<-c(N_age,N_strain_1)
  dim(n_ageoA2)<-c(N_age,N_strain_1)
  dim(n_ageoR2)<-c(N_age,N_strain_1)
  dim(n_ageoE3)<-c(N_age,N_strain_1)
  dim(n_ageoI3)<-c(N_age,N_strain_1)
  dim(n_ageoA3)<-c(N_age,N_strain_1)
  dim(n_ageoR3)<-c(N_age,N_strain_1)
  dim(n_ageoE4)<-c(N_age,N_strain_1)
  dim(n_ageoI4)<-c(N_age,N_strain_1)
  dim(n_ageoA4)<-c(N_age,N_strain_1)
  dim(n_ageoR4)<-c(N_age,N_strain_1)
  dim(n_ageoE12)<-c(N_age,N_strain_2)
  dim(n_ageoI12)<-c(N_age,N_strain_2)
  dim(n_ageoA12)<-c(N_age,N_strain_2)
  dim(n_ageoR12)<-c(N_age,N_strain_2)
  dim(n_ageoE13)<-c(N_age,N_strain_2)
  dim(n_ageoI13)<-c(N_age,N_strain_2)
  dim(n_ageoA13)<-c(N_age,N_strain_2)
  dim(n_ageoR13)<-c(N_age,N_strain_2)
  dim(n_ageoE14)<-c(N_age,N_strain_2)
  dim(n_ageoI14)<-c(N_age,N_strain_2)
  dim(n_ageoA14)<-c(N_age,N_strain_2)
  dim(n_ageoR14)<-c(N_age,N_strain_2)
  dim(n_ageoE23)<-c(N_age,N_strain_2)
  dim(n_ageoI23)<-c(N_age,N_strain_2)
  dim(n_ageoA23)<-c(N_age,N_strain_2)
  dim(n_ageoR23)<-c(N_age,N_strain_2)
  dim(n_ageoE24)<-c(N_age,N_strain_2)
  dim(n_ageoI24)<-c(N_age,N_strain_2)
  dim(n_ageoA24)<-c(N_age,N_strain_2)
  dim(n_ageoR24)<-c(N_age,N_strain_2)
  dim(n_ageoE34)<-c(N_age,N_strain_2)
  dim(n_ageoI34)<-c(N_age,N_strain_2)
  dim(n_ageoA34)<-c(N_age,N_strain_2)
  dim(n_ageoR34)<-c(N_age,N_strain_2)
  dim(n_ageoE4rd)<-c(N_age,N_strain)
  dim(n_ageoI4rd)<-c(N_age,N_strain)
  dim(n_ageoA4rd)<-c(N_age,N_strain)
  dim(n_ageoR4rd)<-c(N_age,N_strain)
  
  
  
  dim(n_muG)<-N_age 
  dim(n_muM)<-N_age
  dim(n_muS)<-N_age
  dim(n_muE)<-c(N_age,N_strain)
  dim(n_muI)<-c(N_age,N_strain)
  dim(n_muA)<-c(N_age,N_strain)
  dim(n_muR)<-c(N_age,N_strain)
  dim(n_muE1)<-c(N_age,N_strain_1)
  dim(n_muI1)<-c(N_age,N_strain_1)
  dim(n_muA1)<-c(N_age,N_strain_1)
  dim(n_muR1)<-c(N_age,N_strain_1)
  dim(n_muE2)<-c(N_age,N_strain_1)
  dim(n_muI2)<-c(N_age,N_strain_1)
  dim(n_muA2)<-c(N_age,N_strain_1)
  dim(n_muR2)<-c(N_age,N_strain_1)
  dim(n_muE3)<-c(N_age,N_strain_1)
  dim(n_muI3)<-c(N_age,N_strain_1)
  dim(n_muA3)<-c(N_age,N_strain_1)
  dim(n_muR3)<-c(N_age,N_strain_1)
  dim(n_muE4)<-c(N_age,N_strain_1)
  dim(n_muI4)<-c(N_age,N_strain_1)
  dim(n_muA4)<-c(N_age,N_strain_1)
  dim(n_muR4)<-c(N_age,N_strain_1)
  dim(n_muE12)<-c(N_age,N_strain_2)
  dim(n_muI12)<-c(N_age,N_strain_2)
  dim(n_muA12)<-c(N_age,N_strain_2)
  dim(n_muR12)<-c(N_age,N_strain_2)
  dim(n_muE13)<-c(N_age,N_strain_2)
  dim(n_muI13)<-c(N_age,N_strain_2)
  dim(n_muA13)<-c(N_age,N_strain_2)
  dim(n_muR13)<-c(N_age,N_strain_2)
  dim(n_muE14)<-c(N_age,N_strain_2)
  dim(n_muI14)<-c(N_age,N_strain_2)
  dim(n_muA14)<-c(N_age,N_strain_2)
  dim(n_muR14)<-c(N_age,N_strain_2)
  dim(n_muE23)<-c(N_age,N_strain_2)
  dim(n_muI23)<-c(N_age,N_strain_2)
  dim(n_muA23)<-c(N_age,N_strain_2)
  dim(n_muR23)<-c(N_age,N_strain_2)
  dim(n_muE24)<-c(N_age,N_strain_2)
  dim(n_muI24)<-c(N_age,N_strain_2)
  dim(n_muA24)<-c(N_age,N_strain_2)
  dim(n_muR24)<-c(N_age,N_strain_2)
  dim(n_muE34)<-c(N_age,N_strain_2)
  dim(n_muI34)<-c(N_age,N_strain_2)
  dim(n_muA34)<-c(N_age,N_strain_2)
  dim(n_muR34)<-c(N_age,N_strain_2)
  dim(n_muE4rd)<-c(N_age,N_strain)
  dim(n_muI4rd)<-c(N_age,N_strain)
  dim(n_muA4rd)<-c(N_age,N_strain)
  dim(n_muR4rd)<-c(N_age,N_strain)
  
  
  #dim(cumu_inc) <- N_age
  dim(seroprev_gii4) <- 7
  dim(inc_day_gii4)<- 5
  dim(inc_day_gii)<- 5
  dim(inc_day_gi3)<- 5
  dim(inc_day_gi)<- 5
  
  dim(inc_day_all_gii4)<- N_age
  dim(inc_day_all_gii)<- N_age
  dim(inc_day_all_gi3)<- N_age
  dim(inc_day_all_gi)<- N_age
  
  dim(inc_day_all)<- 5
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
  dim(foi_gi3) <- N_age
  dim(foi_gi)  <- N_age
  dim(foi_gii4)<- N_age
  dim(foi_gii) <- N_age
  dim(reported_wk)<- N_age
  dim(death_reported_wk_all)<- N_age
  
  
  
  dim(	n_MS	)<-  N_age
    dim(n_bG)   <- N_age
  dim(	n_Mvacc	)<-  N_age
  
  
  dim(n_bM)   <- N_age
  dim(n_bS)   <- N_age
  dim(p_mu)   <- N_age
  dim(p_cfr)   <- N_age
  dim(p_aging)<- N_age
  dim(p_MS)   <- N_age
 dim(p_SE)   <- c(N_age,N_strain)
  dim(p_vacc) <- N_age
  
  dim(pop)<-N_age
  dim(mu) <-N_age
  dim(cfr)<-N_age
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




