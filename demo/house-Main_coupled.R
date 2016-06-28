#
# Main code for partially coupled ABC for household epidemics
#
#

tstartC=proc.time()

#
# Set dataset: SeattleA, SeattleB, Tecumseh1, Tecumseh2
#

data(SeattleA)
DATA=SeattleA

#
# Process data
#

proc = process(DATA)
hsA = proc[[1]]; isA = proc[[2]]; colA = proc[[3]]; rowA = proc[[4]]; xA = proc[[5]]; N = proc[[6]]

#
# Set precision thresholds (epss), infectious period (k) and
# number of iterations (run)
#

epss=c(8,1)
k=0
run=1000

OUTP=PCOUP(DATA,epss,k,run)
OUT=OUTP$OUTPUT
OUTP$simcount

#
# Compute posterior means and standard deviations for partially coupled ABC
#

wei=0
moMG=rep(0,2)
moML=rep(0,2)

Count=length(OUT[,1])
# Number of stored values - note there might be one for each accepted simulation

for(i in 1:Count)
{
wei=wei+Mexp(0,1,0,OUT[i,3],OUT[i,4])
for(j in 1:2)
{
moMG[j]=moMG[j]+Mexp(j,1,0,OUT[i,3],OUT[i,4])
moML[j]=moML[j]+Mexp(0,1,0,OUT[i,3],OUT[i,4])*OUT[i,5]^j
}
}

meanG=moMG[1]/wei
sdG=sqrt(moMG[2]/wei-meanG^2)
meanL=moML[1]/wei
sdL=sqrt(moML[2]/wei-meanL^2)

#
# Computes transformed means and standard deviations of
# q_G = exp(-lambdaG * xA/N); q_L = exp(-lambdaL)
#

A1=1+xA/N
A2=1+2*xA/N

WEIq=exp(-OUT[,3])-exp(-OUT[,4])
moqG=sum(exp(-A1*OUT[,3])-exp(-A1*OUT[,4]))/A1
moqG[2]=sum(exp(-A2*OUT[,3])-exp(-A2*OUT[,4]))/A2

moqL=sum(WEIq*exp(-OUT[,5]))
moqL[2]=sum(WEIq*exp(-2*OUT[,5]))

meanqG=moqG[1]/wei
sdqG=sqrt(moqG[2]/wei-meanqG^2)
meanqL=moqL[1]/wei
sdqL=sqrt(moqL[2]/wei-meanqL^2)

#
# Summarise results
#

# Parameter means and sd
c(meanG,sdG,meanL,sdL)

# Transformed parameters means and sd (compare with Clancy and O'Neill (2008)
# and Neal (2012))
c(meanqG,sdqG,meanqL,sdqL)

# Time taken.

tendC=proc.time()
tendC-tstartC
