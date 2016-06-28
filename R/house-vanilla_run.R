#
# Code for running vanilla ABC for households
#

House_van=function(Xdata,epsil,k,run)
{
OUTPUT=matrix(0,ncol=2,nrow=run)
hs=colSums(Xdata)
isB=rowSums(Xdata)

rowB=seq(0,length(hs))
xB=sum(isB*rowB)

simcount=0
j=0
while(j<run)
{
#print(simcount)
simcount=simcount+1
lambda_G=rexp(1,1)
lambda_L=rexp(1,1)
J=House_SEL(hs,lambda_G,lambda_L,k)
if(sum(abs(J-Xdata))<=epsil[1])
{
isJ=rowSums(J)
if(abs(sum(isJ*rowB)-xB)<=epsil[2])
{
j=j+1
OUTPUT[j,]=c(lambda_G,lambda_L)
print(c(j,simcount))
}
}
}
list(OUTPUT=OUTPUT,simcount=simcount)
}