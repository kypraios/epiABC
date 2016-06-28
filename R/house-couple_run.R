#
# Partially coupled ABC algorithm to obtain run accepted values.
#

PCOUP=function(Xdata,epss,k,run)
{

OUTPUT=matrix(0,ncol=5,nrow=((2*epss[2]+1)*run))
simcount=0
count=0
jj=0

while(jj<run)
{
simcount=simcount+1
lambda_L=rexp(1,1)       # Sample lambda_L
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
OUTPUT[count,]=c(W[ii,2:5],lambda_L) 
print(c(jj,count,simcount))
# Stores values from simulation - these include closeness of simulated epidemic 
# to data, range of lambda_G values and lambda_L
}
}
}
}
}
print(simcount)
list(OUTPUT=OUTPUT[1:count,],simcount=simcount)
}


#PCOUP(DATA,c(20,2),0,1000)


