model{

## Priors
Jmax0 ~ dlnorm(4.7,2.7)             ## maximum electron transport rate prior
alpha0~dnorm(0.25,100)             ##quantum yield  (mol electrons/mole photon) prior
vmax0 ~dlnorm(4.6,2.7)             ## maximum rubisco capacity prior

#Jmax ~ dweibull(2.0,260)          ## maximum electron transport rate prior Serbin 2012
#alpha0 ~ dgamma(2.0,22.0)         ## quantum yield prior Serbin 2012
#vmax0 ~ dweibull(1.7,80)          ## maximum rate of carboxylation prior Serbin 2012

r0 ~ dlnorm(0.75,1.56)             ## leaf respiration prior
#r ~ dweibull(2.0,6.0)             ## broad leaf respiration prior for trees
cp0 ~ dlnorm(1.9,2.7)              ## CO2 compensation point prior
tau ~ dgamma(0.1,0.1)
#TPU  tpu~ dlnorm(3,2.8)             ##tpu

## Constants: Bernacchi et al 2001, PC&E, Table 1
R <- 8.3144621 ## gas constant
r.c <- 18.72
r.H <- 46.39
Vc.c <- 26.35
Vc.H <- 65.33
Vo.c <- 22.98
Vo.H <- 60.11
cp.c <- 19.02
cp.H <- 37.83
cp.ref <- 42.75
Kc.c <- 38.05
Kc.H <- 79.43
Kc.ref <- 404.9
Ko.c <- 20.30
Ko.H <- 36.38
Ko.ref <- 278.4

## Constants: June et al 2004, Funct Plant Bio
Omega <- 18
To <- 35    ## Representative value, would benifit from spp calibration!

## Vcmax BETAs

#RLEAF.V  tau.Vleaf~dgamma(0.1,0.1)          ## add random leaf effects
#RLEAF.V  for(i in 1:nrep){                  
#RLEAF.V   Vleaf[i]~dnorm(0,tau.Vleaf)
#RLEAF.V  }

## alpha BETAs

#RLEAF.A  tau.Aleaf~dgamma(0.1,0.1)
#RLEAF.A  for(i in 1:nrep){                  
#RLEAF.A   Aleaf[i]~dnorm(0,tau.Aleaf)
#RLEAF.A  }

for(i in 1:n) {

r[i]  <- r0 ##B01* exp(r.c - r.H/R/T[i])
cp[i] <- cp0 ##B01* exp(cp.c - cp.H/R/T[i])/cp.ref
Kc.T[i] <- Kc ##B01* exp(Kc.c - Kc.H/R/T[i])/Kc.ref
Ko.T[i] <- Ko ##B01* exp(Ko.c - Ko.H/R/T[i])/Ko.ref
Jmax[i] <- Jmax0 ##J04 * exp(-(T[i]-To)*(T[i]-To)/(Omega*Omega))

alpha[i] <- alpha0 #AFORMULA
al[i]<-(alpha[i]*q[i]/(sqrt(1+(alpha[i]*alpha[i]*q[i]*q[i])/(Jmax[i]*Jmax[i]))))*(pi[i]-cp[i])/(4*pi[i]+8*cp[i])    ## electron transport limited without covariates

vmax.refT[i] <- vmax0 #VFORMULA
vmax[i] <- vmax.refT[i] ##B01* exp(Vc.c - Vc.H/R/T[i])
ae[i]<- vmax[i]*(pi[i]-cp[i])/(pi[i]+Kc.T[i]*(1+po/Ko.T[i]))                                                    ## maximum rubisco limited without covariates

#TPU    ap[i]<-3*tpu                      ## phosphate limited

pmean[i]<-min(al[i], ae[i]) - r[i]      ## predicted net photosynthesis
an[i]~dnorm(pmean[i],tau)            ## likelihood
pA[i] ~ dnorm(pmean[i],tau)          ## prediction
}

foo <- rep[1] + nrep + T[1]                ## prevent warnings
}
