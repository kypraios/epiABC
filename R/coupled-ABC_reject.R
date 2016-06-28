#
# ABC rejection sampler for a fixed number of accepted values.
#

ABCrej=function(N,mstar,epsil,run)
{
samp=rep(0,run)
simcount=0
i=0

while(i<run)
{
simcount=simcount+1
lambda=rexp(1)
m=SIR_sim(N,lambda,1)
if(abs(m-mstar)<=epsil)
{
i=i+1
samp[i]=lambda
#print(c(i,simcount))
}
}
print(simcount) # Number of simulations required to obtain the sample.
list(samp=samp,simcount=simcount)
}