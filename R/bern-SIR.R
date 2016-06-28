bernSIR=function(n,beta,gamma,p)
{
t=0
MAT=matrix(rbinom(n^2,1,p),ncol=n)
# Matrix of edges - non-symmetric but ok because we will only use MAT[i,j] if
# i infected before j or MAT[j,i] if j infected before i.
rowM=rowSums(MAT)
I=rep(0,n)
I[1]=1
# Set individual 1 infectious and everybody else susceptible.
output=rep(0,n) # Recovery times
count=0  # Number of recoveries observed
while(sum(I==1)>0)
{
rec=sum(I==1)
infe=sum(rowM[I==1])
t=t+rexp(1,(gamma*rec+beta*infe))
u=runif(1)
if(u<=(beta*infe/(gamma*rec+beta*infe)))
{
S=rep(0,n)
S[I==1]=rowM[I==1]
K=sample(n,1,replace=TRUE,prob=S)
J=sample(n,1,replace=TRUE,prob=MAT[K,])
if(I[J]==0) I[J]=1
}
if(u>(beta*infe/(gamma*rec+beta*infe)))
{
S=rep(0,n)
S[I==1]=1
K=sample(n,1,replace=TRUE,prob=S)
I[K]=2
count=count+1
output[count]=t
}
}
#output[1:count]
list(output=output,count=count)
}
