#
# Code for simulating a household epidemic using the Sellke construction
#

House_SEL=function(hs,lambda_G,lambda_L,k)
{
HH=length(hs)
ks=seq(1,length(hs))

n=rep(ks,hs) # size of each household

m=length(n) # Number of households
N=sum(n) # Population size
NS=N

sev=0
threshold=0

ni=rep(0,length(n)) # infectives
ns=n #susceptibles

R=rexp(N)

while(threshold<=(lambda_G*sev))
{
kk=sample(m,1,replace=T,ns)

hou_epi=House_epi(ns[kk],k,lambda_L)
ns[kk]=ns[kk]-hou_epi[1]
sev=sev+hou_epi[2]

ni[kk]=n[kk]-ns[kk]
#NS=NS-hou_epi[1]

NS=sum(ns)

if(NS>0) threshold=threshold+rexp(1,(NS/N))
if(NS==0) threshold=2*lambda_G*sev
}
OUT=matrix(0,ncol=HH,nrow=(HH+1))
for(i in 1:(HH+1))
{
for(j in 1:HH)
{
OUT[i,j]=sum(n[ni==(i-1)]==j)
}
}
OUT
}
