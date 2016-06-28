#
# Code for running SMC-ABC for partially coupled household epidemics.
#

PCOUP_SMC=function(Xdata,epss,k,run,OX)
{

OUTPUT=matrix(0,ncol=6,nrow=((2*epss[2]+1)*run))
simcount=0
count=0
jj=0

#
# Setting standard deviation for importance sampling
#

VL=sum(OX[,3]^2*OX[,4])-sum(OX[,3]*OX[,4])^2
sL=sqrt(2*VL)
print(sL)

#
# Number of samples stored from run acceptances.
#

cox=length(OX[,1])

while(jj<run)
{
simcount=simcount+1
LA=sample(cox,1,replace=TRUE,OX[,4])
lambda_L=rnorm(1,OX[LA,3],sL)       # Sample lambda_L
if(lambda_L>0)
{
J=House_COUP(Xdata,epss[2],lambda_L,k) # Run coupled simulations
W=J[J[,4]<J[,5],]        # W contains successful simulations (infect close to xA individuals)
if(length(W)==5) W=matrix(W,ncol=5,byrow=T)
if(length(W[,1])>0)
{
if(min(W[,2])<=epss[1])
{
jj=jj+1
for(ii in 1:length(W[,1]))
{
if(W[ii,2]<=epss[1])
{
count=count+1
OUTPUT[count,]=c(W[ii,2:5],lambda_L,(exp(-lambda_L)/sum(dnorm(lambda_L,OX[,3],sL)))) 
#print(c(jj,count,simcount))
# Stores values from simulation - these include closeness of simulated epidemic 
# to data, range of lambda_G values and lambda_L
}
}
}
}
}
}
print(simcount)
list(OUTPUT=OUTPUT[1:count,],simcount=simcount)
}

