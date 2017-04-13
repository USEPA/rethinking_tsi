model{
  for (i in 1:1078){
    
    for (num_cutpt in 1:3){
      Z[i,num_cutpt] <- ( alpha_SD*SD[i] 
#                           + alpha_N*mu_N[i]
#                           + alpha_P*mu_P[i]
                          + alpha_N*Nitrogen[i] 
                          + alpha_P*Phosphorus[i]
                          + alpha_E*Eleveation[i]
                          - C[num_cutpt])
      # / s
      
      Q[i,num_cutpt] <- 1/(1+exp(-Z[i,num_cutpt]))
    }
    
    P[i,1] <- max(min(1 - Q[i,1],1),0)
    P[i,2] <- Q[i,1] - Q[i,2]
    P[i,3] <- Q[i,2] - Q[i,3]
    P[i,4] <- max(min(Q[i,3],1),0)
    
    TS[i] ~ dcat(P[i,])
    
#     Phosphorus[i] ~ dnorm(mu_P[i], tau_P)
#     mu_P[i] <- gamma_L*Latitude[i]
#     +  gamma_E*Evergreen[i]
#     +  gamma_ER[Eco_Region[i]]
#     
#     Nitrogen[i] ~ dnorm(mu_N[i], tau_N)
#     mu_N[i] <- beta_L*Latitude[i]
#     +  beta_E*Evergreen[i]
#     +  beta_ER[Eco_Region[i]]
  }
  
  # PRIORS
  # POLR Coefficients
  alpha_SD ~ dnorm(0,.0001)
  alpha_N ~ dnorm(0,.0001)
  alpha_P ~ dnorm(0,.0001)
  alpha_E ~ dnorm(0,.0001)
  # Nitrogen Multilevel Coefficients
#   tau_N <- pow(sigma_N_hat, -2)
#   sigma_N_hat ~ dunif(1,100)
#   beta_L ~ dnorm(0,.0001)
#   beta_E ~ dnorm(0,.0001)
#   for (j in 1:9){
#     beta_ER[j] ~ dnorm(0,.0001)
#   }
#   
  # Phosphorus Multilevel Coefficients
#   tau_P <- pow(sigma_P_hat, -2)
#   sigma_P_hat ~ dunif(1,100)
#   gamma_L ~ dnorm(0,.0001)
#   gamma_E ~ dnorm(0,.0001)
#   for (j in 1:9){
#     gamma_ER[j] ~ dnorm(0,.0001)
#   }
  # s ~ dlnorm(mu.log.s, tau.log.s)
  # mu.log.s ~ dnorm(0,0.0001)
  # tau.log.s <- pow(sigma.log.s,-2)
  # sigma.log.s ~ dunif(0,1000)
  
  for (i in 1:3) {
    cutpt_raw[i] ~ dnorm(0, .0001)
  }
  C <- sort(cutpt_raw)
  
  # C[1] ~ dnorm(0,.0001)I(,C[2])
  # C[2] ~ dnorm(0,.0001)I(C[1],C[3])
  # C[3] ~ dnorm(0,.0001)I(C[2],)
  
}