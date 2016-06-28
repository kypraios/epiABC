#
# Computing mean and variance
#

dpow=function(MM,k)
{
(sum(MM[,2]**k)-sum(MM[,1]**k))/k
}

meanC=function(MM)
{
dpow(MM,2)/dpow(MM,1)
}

varC=function(MM)
{
dpow(MM,3)/dpow(MM,1)-meanC(MM)**2
}

sdC=function(MM)
{
varC(MM)**(1/2)
}

