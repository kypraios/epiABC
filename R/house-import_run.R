#
# Code for running ABC for households with importance sampling
#

House_imp=function(Xdata,epsil,k,run,OX)
{
OUTPUT=matrix(0,ncol=3,nrow=run)
hs=colSums(Xdata)
isB=rowSums(Xdata)

rowB=seq(0,length(hs))
xB=sum(isB*rowB)

simcount=0
j=0

#
# Compute variances for the pertubations. Bivariate Gaussian
# 

meanG=sum(OX[,1]*OX[,3])/sum(OX[,3])
meanL=sum(OX[,2]*OX[,3])/sum(OX[,3])
vG=sum(OX[,1]^2*OX[,3])/sum(OX[,3])-meanG^2
vL=sum(OX[,2]^2*OX[,3])/sum(OX[,3])-meanL^2
vLG=sum(OX[,1]*OX[,2]*OX[,3])/sum(OX[,3])-meanG*meanL
vaR=2*matrix(c(vG,vLG,vLG,vL),ncol=2)
Sinv=solve(vaR)
sG=sqrt(2*vG)
sLL=sqrt(2*(vL-vLG^2/vG))
sAL=vLG/vG

while(j<run)
{
simcount=simcount+1
LA=sample(run,1,replace=TRUE,OX[,3])
lambda_G=rnorm(1,OX[LA,1],sG)
lambda_L=rnorm(1,(OX[LA,2]+sAL*(lambda_G-OX[LA,1])),sLL)

if((lambda_G>0)&(lambda_L>0))
{
J=House_SEL(hs,lambda_G,lambda_L,k)
if(sum(abs(J-Xdata))<=epsil[1])
{
isJ=rowSums(J)
if(abs(sum(isJ*rowB)-xB)<=epsil[2])
{
j=j+1
weiZ=0
for(i in 1:run)
{
xDIF=c(lambda_G,lambda_L)-OX[LA,1:2]
mult=t(xDIF)%*%Sinv%*%xDIF
weiZ=weiZ+exp(-mult[1,1]/2)
}
weiG=exp(-(lambda_L+lambda_G))/weiZ
OUTPUT[j,]=c(lambda_G,lambda_L,weiG)
print(c(j,simcount))
}
}
}
}
list(OUTPUT=OUTPUT,simcount=simcount)
}