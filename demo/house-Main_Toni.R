#
# ABC code for running ABC using Toni et al. SMC
#

tstartC=proc.time()

#
# Set dataset: SeattleA, SeattleB, Tecumseh1, Tecumseh2
#

data(Tecumseh1)
DATA=Tecumseh1

#
# Process data
#

proc = process(DATA)
hsA = proc[[1]]; isA = proc[[2]]; colA = proc[[3]]; rowA = proc[[4]]; xA = proc[[5]]; N = proc[[6]]

#
# Code for running coupled algorithm and selecting which simulations to keep.
#

#
# Set thresholds
#

Eps1=c(100,70,50)
Eps2=c(10,6,4)

#Eps1=c(20,12,8)
#Eps2=c(3,2,1)

TT=length(Eps1)
EpsM=matrix(c(Eps1,Eps2),ncol=2,byrow=F)

# Set infectious period (k) and  number of iterations (run)
k=0
run=1000

#
# Initial run
#

epss=EpsM[1,]

OUTP=House_van(DATA,EpsM[1,],k,run)
OUT=OUTP$OUTPUT
simT=OUTP$simcount

OL=matrix(0,ncol=3,nrow=run)
OL[,1:2]=OUT
OL[,3]=rep(1,run)

#
# Subsequent runs using importance sampling
#

if(TT>1)
{
for(t in 2:TT)
{
simT
OUTP=House_imp(DATA,EpsM[t,],k,run,OL)
OL=OUTP$OUTPUT
simT=simT+OUTP$simcount
}
}

OUT=OL

#
# Moment calculations for SMC-ABC output
#

weiA=sum(OUT[,3])

meanG=sum(OUT[,1]*OUT[,3])/weiA
sdG=sqrt(sum(OUT[,1]^2*OUT[,3])/weiA-meanG^2)
meanL=sum(OUT[,2]*OUT[,3])/weiA
sdL=sqrt(sum(OUT[,2]^2*OUT[,3])/weiA-meanL^2)

#
# Moment calculations for transformed variables
#

meanqG=sum(exp(-OUT[,1]*xA/N)*OUT[,3])/weiA
sdqG=sqrt(sum(exp(-2*OUT[,1]*xA/N)*OUT[,3])/weiA-meanqG^2)
meanqL=sum(exp(-OUT[,2])*OUT[,3])/weiA
sdqL=sqrt(sum(exp(-2*OUT[,2])*OUT[,3])/weiA-meanqL^2)

#
# Summarise results
#

c(meanG,sdG,meanL,sdL)

# Transformed parameters.
c(meanqG,sdqG,meanqL,sdqL)

simT

tendC=proc.time()
tendC-tstartC
