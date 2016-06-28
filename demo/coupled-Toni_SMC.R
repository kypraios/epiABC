#
# Code for implementing the Toni et al. alogrithm with predefined levels of
# approximation.
#
# Implementation for the Abakiliki data

tstartC=proc.time()

#
#

# Data

N=120
mstar=30

# Predefined thresholds.

epsil=c(10,5,2,1,0)
epsil=c(10,0)
TT=length(epsil)

# Estimate of posterior mean and variance at each threshold.
EstMean=rep(0,TT)
EstSD=rep(0,TT)
# Total number of simulations
SIMtotal=rep(0,TT)

# Number of particles.

run=10000

output=ABCrej(N,mstar,epsil[1],run)
samp=output$samp
simTotal=output$simcount

EstMean[1]=mean(samp)
EstSD[1]=sd(samp)
SIMtotal[1]=simTotal

# Initial weights all equal
weight=rep(1,run)

for(t in 2:TT)
{
output=ABCimp(N,mstar,epsil[t],run,samp,weight)
samp=output$samp
weight=output$weight
simTotal=simTotal+output$simcount

SIMtotal[t]=simTotal
EstMean[t]=sum(samp*weight)/sum(weight)
EstSD[t]=sqrt(sum(samp^2*weight)/sum(weight)-EstMean[t]^2)
}

#
#

EstMean
EstSD
SIMtotal


tendC=proc.time()
tendC-tstartC
