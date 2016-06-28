# vanilla ABC for the Abakaliki Data

#Day   0 1 2 3 4 5 6 7
#Cases 1 0 4 2 3 3 10 5

gastro.removal.times <- c(0, 2,2,2,2,3,3,4,4,4,5,5,5,rep(6, 10), rep(7, 5))
breaks.gastro <- c(0, 1, 2, 3, 4, 5, 6, 7, Inf)
gastro.binned.data <- hist(gastro.removal.times, breaks=breaks.gastro, plot=FALSE)$count
gastro.binned.data
# c(1,0,4,2,3,3,10,5)

ptm <- proc.time()
out <- abcSIR.binned(obs.data = gastro.binned.data, breaks.data = breaks.gastro, obs.duration = max(gastro.removal.times), N = 89, epsilon = 10, prior.param = c(0.1,0.1), samples = 500)
atm <- proc.time() - ptm
cat("It took", atm[3]/60, "mins")

hist(out[,1]);
hist(out[,2]);

hist(out[,1]/out[,2])
mean(out[,1]/out[,2])

mean(out[,1])
sd(out[,1])

mean(out[,2])
sd(out[,2])




#####
##[1] 500
##[1] 120667317
