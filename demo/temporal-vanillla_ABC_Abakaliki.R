# vanilla ABC for the Abakaliki Data

abakaliki.removal.times <- cumsum(c(0,13,7,2,3,0,0,1,4,5,3,2,0,2,0,5,3,1,4,0,1,1,1,2,0,1,5,0,5,5))
breaks.abakaliki <- c(0, 13, 26, 39, 52, 65, 78, Inf)
abakaliki.binned.data <- hist(abakaliki.removal.times, breaks=breaks.abakaliki, plot=FALSE)$count

ptm <- proc.time()
out <- abcSIR.binned(obs.data = abakaliki.binned.data, breaks.data = breaks.abakaliki, obs.duration = max(abakaliki.removal.times), N = 120, epsilon = 11, prior.param = c(0.1,0.1), samples = 500)
atm <- proc.time() - ptm
cat("It took", atm[3]/60, "mins")

#write.table(out, file="vanilla_ABC_epsilon_10_prior_0.1_0.1.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(out, file="vanilla_ABC_epsilon_11_prior_0.1_0.1.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)

# values to compare
# posterior mean for beta using a GIHM: 9.5*(10^(-4))
# posterior variance for beta using a GIHM: 6.6*(10^(-8))
# posterior mean for gamma using a GIHM 0.10 
# posterior variance for gamma using GIHM: 8.3*10^(-4)


#hist(out[,1]); abline(v=120*9.5*(10^(-4)),col=2)
#hist(out[,2]); abline(v=0.1, col=2)

mean(out[,1])
sd(out[,1])

mean(out[,2])
sd(out[,2])
