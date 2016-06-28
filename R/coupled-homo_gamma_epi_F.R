#
# Homogeneously mixing SIR code
#

# 
# Simulates epidemics with a Gamma(k,k) infectious period.
# output gives the set of lambda parameters consistent with the data (m out of n infected)
# Only successful simulations kept 
# k=0 constant infectious period

epigammaF=function(n,m,run,k)
{
output=matrix(0,ncol=2,nrow=run)
count=0
j=0
while(j<run)
{
count=count+1
t=thres(n)
if(k>0) q=rgamma(n,k,k)
if(k==0) q=rep(1,n)
y=0
for(i in 1:(n-1)) y[i]=t[i]/sum(q[1:i])
q=max(y[1:(m-1)])
if(q<y[m])
{
j=j+1
output[j,]=c(q,y[m])
}
}
print(count)
output
}

