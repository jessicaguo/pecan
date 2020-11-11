model{
  
  ## Priors
  #alpha0 ~ dlnorm(-3.21,3.7) 	    	 ## initial slope of photosynthesis light response
  alpha0 ~ dgamma(90.9, 1580) 
  vmax0~dlnorm(3,3)                  ## maximum rubisco capacity
 # r0 ~ dlnorm(-0.2,2.8)              ## leaf respiration
  k0 ~ dlnorm(11.5, 2.8)             ## initial slope of photosynthetic CO2 response
  tau ~ dgamma(0.1,0.1)
  r0 ~ dbeta(2,2)
  # r0 ~ dnorm(rd, tau.r0)
  # tau.r0 ~ dgamma(0.1, 0.1)
  
  ## alpha BETAs
  
  #RLEAF.A  tau.Aleaf~dgamma(0.01,0.01)
  #RLEAF.A  for(i in 1:nrep){                  
  #RLEAF.A   Aleaf[i]~dnorm(0,tau.Aleaf)
  #RLEAF.A   A[i] <- Aleaf[i] + alpha0
  #RLEAF.A  }
  
  ## Vcmax BETAs
  
  #RLEAF.V  tau.Vleaf~dgamma(0.01,0.01)          
  #RLEAF.V  for(i in 1:nrep){                  
  #RLEAF.V   Vleaf[i]~dnorm(0,tau.Vleaf)
  #RLEAF.V   V[i] <- Vleaf[i] + vmax0
  #RLEAF.V  }
  
  ## k BETAs
  
  #RLEAF.K  tau.Kleaf~dgamma(0.01,0.01)
  #RLEAF.K  for(i in 1:nrep){                  
  #RLEAF.K   Kleaf[i]~dnorm(0,tau.Kleaf)
  #RLEAF.K   K[i] <- Kleaf[i] + k0
  #RLEAF.K  }
  
  ## r BETAs
  
  #RLEAF.R  tau.Rleaf~dgamma(0.01,0.01)
  #RLEAF.R  for(i in 1:nrep){                  
  #RLEAF.R   Rleaf[i]~dnorm(0,tau.Rleaf)
  #RLEAF.R   R[i] <- Rleaf[i] + r0*rd
  #RLEAF.R  }
  
  for(i in 1:n){ 
    
    ## light limited
    alpha[i] <- alpha0 #AFORMULA
    al[i]<- alpha[i]*q[i]  
    
    ## CO2 limited
    k[i] <- k0 #KFORMULA
    ac[i]<-k[i]*pi[i]/100000     
    
    ## Rubisco limited
    vmax.refT[i] <- vmax0 #VFORMULA
    vmax[i] <- vmax.refT[i] ##B01* exp(Vc.c - Vc.H/R/T[i])
    ae[i] <- vmax[i]    
    
    ## Dark respiration
    r[i] <- r0*rd #RFORMULA
    
    pmean[i]<-min(min(al[i],ac[i]),ae[i])-r[i]      ## predicted net photosynthesis
    an[i] ~ dnorm(pmean[i],tau)                  ## likelihood
    pA[i] ~ dnorm(pmean[i],tau)                  ## prediction
  }
}