#
#  ABC code for running vanilla ABC
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

#
# Run vanilla ABC
#

OUTP=House_van(DATA,epss,k,run)
OUT=OUTP$OUTPUT

#
# Moment calculations
#

meanG=mean(OUT[,1])
sdG=sd(OUT[,1])
meanL=mean(OUT[,2])
sdL=sd(OUT[,2])

#
# Transformed moments
#

meanqG=mean(exp(-OUT[,1]*xA/N))
sdqG=sd(exp(-OUT[,1]*xA/N))
meanqL=mean(exp(-OUT[,2]))
sdqL=sd(exp(-OUT[,2]))

#
# Summarise results
#

c(meanG,sdG,meanL,sdL)

# Transformed parameters.
c(meanqG,sdqG,meanqL,sdqL)

OUTP$simcount

tendC=proc.time()
tendC-tstartC
