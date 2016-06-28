#
# Homogeneously mixing SIR code
#

#
# Computes successive thresholds
#

thres=function(n)
{
thres=rep(0,(n-1))
thres[1]=rexp(1,(n-1)/n)
for(i in 2:(n-1)) thres[i]=thres[i-1]+rexp(1,(n-i)/n)
thres
}
