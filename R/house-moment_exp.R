#
# Moment calculator for Exp(\lambda) prior.
#

Mexp=function(k,lambda,U,a,b)
{
if(k==0) U=U+exp(-lambda*a)-exp(-lambda*b)
if(k>0) U=U+(a^k*exp(-lambda*a)-b^k*exp(-lambda*b))+(k/lambda)*Mexp((k-1),lambda,U,a,b)
U
}

#Mexp(3,2,0,1,3)


