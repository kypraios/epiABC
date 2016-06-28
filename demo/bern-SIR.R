sin=proc.time()

N=1000
OUT=matrix(0,ncol=89,nrow=N)
PARA=matrix(0,ncol=5,nrow=N)

for(i in 1:N)
{
p=runif(1)
beta=rexp(1,2)
gamma=rgamma(1,2,1)
#print(PARA[i,])
Xout=bernSIR(89,beta,gamma,p)
OUT[i,]=Xout$output
PARA[i,]=c(p,beta,gamma,Xout$count,OUT[i,Xout$count])
}

sout=proc.time()
sout-sin
