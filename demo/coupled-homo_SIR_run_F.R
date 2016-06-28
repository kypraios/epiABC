#
# Homogeneously mixing SIR code
#

tstartC=proc.time()

#
# Code for Abikiliki data: n=120, m=30
#

run=10000

# Gamma distributed infectious period
wk=epigammaF(120,30,run,1)

#
#meanC(wk)
#sdC(wk)

moM=rep(0,3)

for(i in 1:run)
{
for(j in 1:3)
{
moM[j]=moM[j]+Mexp((j-1),1,0,wk[i,1],wk[i,2])
}
}

moM[2]/moM[1]
sqrt(moM[3]/moM[1]-(moM[2]/moM[1])^2)


tendC=proc.time()
tendC-tstartC
