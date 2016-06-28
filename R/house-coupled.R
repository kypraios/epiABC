#  Partially coupled ABC algorithm for household epidemics
#  lambda_L is drawn from the prior (or however)
#  Code finds lambda_G values consistent with the data
#
#  Input: Xdata - Epidemic data to compare simulations with.
#  epsil - Max distance between simulated and observed final size for a simulation 
#  to be accepted. (Tighter control on distance after simulations straightforward).
#  lambda_L - local infection (household) rate
#  k - Gamma(k,k) infectious period with k=0 a constant infectious period.
#

House_COUP=function(Xdata,epsil,lambda_L,k)
{

hsA=colSums(Xdata)   # hsA[i] - Number of households of size i
isA=rowSums(Xdata)   # isA[i] - Number of households with i-1 infectives
colA=seq(1,length(hsA))
rowA=seq(0,length(hsA))

xA=sum(isA*rowA)   # Final size

HH=length(hsA) # HH maximum household size
ks=seq(1,HH) 

n=rep(ks,hsA) # household data

m=length(n) # Number of households
N=sum(n) # Population size
NS=N     # Number of susceptibles

sev=0       # Running tally of severity (sum of infectious periods)
threshold=0 # Running tally of (global) threshold required for the next infection

ni=rep(0,length(n)) # infectives (per household)
ns=n                # susceptibles (per household)

OUT=matrix(0,ncol=HH,nrow=(HH+1))
OUT[1,]=hsA        # Epidemic data in the same form as Xdata
                   # Start with everybody susceptible

DISS=matrix(0,ncol=5,nrow=(2*epsil+1)) 
# Matrix for collecting epidemics infecting within epsil of xA infectives.

SEVI=matrix(0,ncol=3,nrow=N)
# Matrix to keep track of number of infectives, severity and threshold.

ys=0    # number of infectives
count=0 # number of global infections taking place.
        # first global infection is the introductory case. 

while(ys<=(xA+epsil)) 
# Only need to consider the epidemic until xA+epsil infections occur.
# We simulate successive global infections (should they occur) with associated
# local (within household epidemics)
# For the count+1 global infection to take place, we require that 
# for k=1,2,..., count;  lambda_G * severity from first k infectives is larger
# than the k^th threshold
{
count=count+1
kk=sample(m,1,replace=T,ns) # Choose household for next global infection
# The first kk chooses the initial infective.

OUT[(ni[kk]+1),n[kk]]=OUT[(ni[kk]+1),n[kk]]-1

hou_epi=House_epi(ns[kk],k,lambda_L)
# Simulate a household epidemic among the remaining susceptibles in the household

ns[kk]=ns[kk]-hou_epi[1]
ni[kk]=n[kk]-ns[kk]
# Update household kk data (susceptibles and infectives)

OUT[(ni[kk]+1),n[kk]]=OUT[(ni[kk]+1),n[kk]]+1
# Update the state of the population following the global infection and 
# resulting household epidemic


NS=sum(ns)
threshold=threshold+rexp(1,(NS/N)) 
# Generates a threshold for the next global infection

ys=ys+hou_epi[1] # Update number infected
sev=sev+hou_epi[2] # Update severity
SEVI[count,]=c(ys,sev,threshold)

if(abs(ys-xA)<=epsil)
# If the number infected is close to xA, we check what value of lambda_G 
# would be needed for an epidemic of the desired size. 
# Note that in many cases no value of lambda_G will result in an epidemic 
# close to xA. 
{
dist=sum(abs(OUT-Xdata))   # Measures distance between simulated data and observed
TT=SEVI[1:count,3]/SEVI[1:count,2] # Ratio of threshold to severity
Tlow=max(TT[1:(count-1)])   #  Tlow is the minimum lambda_G which leads to 
                            #  at least count global infections
Thi=max(TT[1:(count)])      #  Thi is the maximum lambda_G which leads to 
                            #  at most count global infections
# In other words, lambda_G between Tlow and Thi leads to exactly count 
# global infections
DISS[(ys+1-(xA-epsil)),]=c(1,dist,abs(ys-xA),Tlow,Thi)
# DISS stores key information on simulations with ys close to xA.
}
}

DISS     # Output from simulations
         # Key interest in rows where Tlow < Thi (column 4 < column 5)
}

