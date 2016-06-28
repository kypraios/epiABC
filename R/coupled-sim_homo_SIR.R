#
# R code for simulating SIR epidemic
#
# Gamma distributed infectious period  - Gamma (k,k)
# (k=0 - constant infectious period)
#

SIR_sim=function(N,lambda,k)
{
S=N-1
y=1
while(y>0)
{
  y=y-1
  if(k==0) I=1
  if(k>0) I=rgamma(1,k,k)
  Z=rpois(1,(lambda*I))
  if(Z>0)
  {
    for(j in 1:Z)
    {
      u=runif(1)
      if(u<(S/N))
      {
        S=S-1
        y=y+1
      }
    }
  }
}
N-S
}
