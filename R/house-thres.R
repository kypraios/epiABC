# Code for setting (local) thresholds in a household of size n.
# 
# Note local thresholds not required for households of size 1.
# Infection rate does not depend upon households size.
#

thresH=function(n)
{
thres=rep(0,(n-1))
thres[1]=rexp(1,(n-1))
if(n>2) 
{
for(i in 2:(n-1)) thres[i]=thres[i-1]+rexp(1,(n-i))
}
thres
}