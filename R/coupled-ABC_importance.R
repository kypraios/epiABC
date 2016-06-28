#
# ABC importance sampling for a fixed number of accepted values.
# Use this in implementing the Toni et al. algorithm.
#

ABCimp=function(N,mstar,epsil,run,sampold,weightold)
# sampold and weightold are the lambda values and importance weight from the
# previous run.
{
samp=rep(0,run)
weight=rep(0,run)
simcount=0
i=0

V=sum(sampold^2*weightold)/sum(weightold)-sum(sampold*weightold)^2/sum(weightold)^2
ss=sqrt(2*V)
print(ss)
#
# Beaumont et al. (2009) for the choice of ss.


while(i<run)
{
simcount=simcount+1
ll=sample(sampold,1,replace=TRUE,weightold)
lambda=rnorm(1,ll,ss)
if(lambda>0)
{
m=SIR_sim(N,lambda,1)
if(abs(m-mstar)<=epsil)
{
i=i+1
samp[i]=lambda
weight[i]=exp(-lambda)/sum(weightold*dnorm(lambda,sampold,ss))
# Computes the importance weight for the parameter based upon an Exp(1) prior.
#print(c(i,simcount))
}
}
}
print(simcount) # Number of simulations required to obtain the sample.
list(samp=samp,weight=weight,simcount=simcount)
}