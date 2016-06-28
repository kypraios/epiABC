# Code for simulating a household epidemic.
#
# 1 - infective; n-1 - susceptible
# Infectious periods Gamma(k,k) - k=0 implies constant infectious period =1.
# Infection rate lamdba_L
#
# Output: number infected and severity (sum of infectious periods)
#

House_epi=function(n,k,lambda_L)
{

# Household of size 1 no local infection

if(n==1)
{
i=1
if(k==0) sev=1
if(k>0) sev=rgamma(1,k,k)
}

# Household of size n>1 local infection possible

if(n>1)
{
t=thresH(n)    # Local thresholds
if(k==0) q=rep(1,n)    # Constant infectious period
if(k>0) q=rgamma(n,k,k) # Gamma infectious period
t[n]=2*lambda_L*sum(q)  # t[n] only used if everybody in the household infected
                        # Upper bound to ensure t[n]>lambda_L*sum(q)
i=0
test=0
while(test==0)
{
i=i+1
if(t[i]>(lambda_L*sum(q[1:i])))
# Epidemic in the household stops when the infection from the first i infectives 
# is insufficient to infect the i+1^st individual 
{
test=1
sev=sum(q[1:i])
}
}
}
c(i,sev)
# Output number infected and severity
}

