# Simulate a Markovian SIR model
# This function assumes 1 initial infective and N-1 initially susceptibles
# Per-person infection rate is beta/N

simSIR <- function(N, beta, gamma) {
   
   # initial number of infectives and susceptibles;
   I <- 1
   S <- N-1;
   
   # recording time;
   t <- 0;
   times <- c(t);
   
   # a vector which records the type of event (1=infection, 2=removal)
   type <- c(1);
   
   while (I > 0) {
      
      # time to next event;
      t <- t + rexp(1, (beta/N)*I*S + gamma*I);
      times <- append(times, t);
      
      if (runif(1) < beta*S/(beta*S + N*gamma)) {
         # infection
         I <- I+1;
         S <- S-1;
         type <- append(type, 1);
      }
      else {
         #removal
         I <- I-1
         type <- append(type, 2);
      }
   }
   
   # record the removal time, the final size, and the durection of the epidemic. 
   # note that we scale the removal times such that the first removal is a time zero.
   #
   res <- list("removal.times" = times[type==2] - min(times[type==2]) , "final.size" = N-S, "T" = times[length(times)])
   
   res
}

# this is function that simulates an outbreak of a predetermined final size -- to be used for the
# Abakaliki data
simSIR.constrained <- function(N, beta, gamma, final.size) {
  check <- 0
  while (check==0) {
    out <- simSIR(N, beta, gamma)
    if (out$final.size >= final.size) {
      res <- out$removal.times
      check <- 1
    }
  }
  return(res)
}

abcSIR <- function(obs.data, N, epsilon, prior.param, samples) {
   
   # first retrieve the final size of the observed data
   final.size.obs <- length(obs.data)
   
   # matrix to store the posterior samples
   post.samples <- matrix(NA, ncol = 2, nrow=samples)
   
   i = 0
   while (i < samples) {
      
      # draw from the prior distribution
      beta <- rexp(1, prior.param[1])
      gamma <- rexp(1, prior.param[2])
      
      # simulate data
      sim.data <- simSIR(N, beta, gamma)
      
      #check if the final size matches the observedata
      if (sim.data$final.size == final.size.obs) {
         d <- sum((obs.data - sim.data$removal.times)^2)

         if (d < epsilon){
            i <- i + 1
            post.samples[i,] <- c(beta, gamma)
            
         }
         
      }
   }
   post.samples
}


abcSIR.binned <- function(obs.data.binned, breaks.data, obs.duration, N, epsilon, prior.param, samples) {
   
   # first retrieve the final size of the observed data
   final.size.obs <- length(obs.data.binned)
   
   # matrix to store the posterior samples
   post.samples <- matrix(NA, ncol = 2, nrow=samples)
   
   K = 0
   
   i = 0
   while (i < samples) {
      
      # counter
      K = K + 1
      
      # draw from the prior distribution
      beta <- rexp(1, prior.param[1])
      gamma <- rexp(1, prior.param[2])
      
      # simulate data
      sim.data <- simSIR(N, beta, gamma)
      sim.duration <- sim.data$T
      sim.data.binned <- hist(sim.data$removal.times, plot=FALSE, breaks=breaks.data)$counts
      
      #check if the final size matches the observedata
      d <- sqrt( sum((obs.data.binned - sim.data.binned)^2) + ((obs.duration - sim.duration)/50)^2 )
         
      if (d < epsilon){
         i <- i + 1
         print(i)
         post.samples[i,] <- c(beta, gamma)
         
      }
   }
   print(K)
   post.samples
}

simSIR.discrete <- function(N, lambda, gamma, T) {
   
   # change the rate to avoidance probability
   q <- exp(-lambda/N)
   
   # initialisation
   t.vec <- seq(1, T)
   
   It.vec <- rep(NA, length = T)
   St.vec <- rep(NA, length = T)
   Rt.vec <- rep(NA, length = T)
   
   It.vec[1] <- 1
   St.vec[1] <- N - 1
   Rt.vec[1] <- 0
   
   # sample infectious period for the initially infective individual
   inf.per <- rgeom(1, gamma) + 1; 
   
   # Yt.vec keeps track of the number of people being removed on each day
   Yt.vec <- rep(0, length = T)
   
   # assing the value to Yt which corresponds to when the initially infective individual gets removed
   Yt.vec[inf.per + 1] <- Yt.vec[inf.per + 1] + 1  
   
   t <- 1
   while (t < T) {
      
      t <- t + 1
      
      # simulate the number of new infections for next day t
      if (It.vec[t-1] > 0) {
         new.inf <- rbinom(1, St.vec[t-1], prob =  1 - q^It.vec[t-1])
         St.vec[t] <- St.vec[t-1] - new.inf
         It.vec[t] <- It.vec[t-1] + new.inf - Yt.vec[t] 
         # when update the vector It.vec we need to take into account any individuals who will be removed then
         
         # simulate for each individual who got infected their infectious period
         if (new.inf > 0) {
            for (j in 1:new.inf) {
               inf.per <- rgeom(1, gamma) + 1
               loc <- min(t + 0 + inf.per, T)
               Yt.vec[loc] <- Yt.vec[loc] + 1
            }
         }            
      }
      else {
         It.vec[t] <- It.vec[t-1]
         St.vec[t] <- St.vec[t-1]
      }
   }
   
   Rt.vec <- cumsum(Yt.vec)
   res <- rbind(It.vec, St.vec, Rt.vec, Yt.vec)
   if (sum(apply(res[-4,], 2, sum) - rep(N, T))!=0) stop("error")
   
   out <- list("pop"=res, "final.size"=sum(res[4,]))
   
   return(out)
   
}


abcSIR.discrete <- function(obs.data, N, T, epsilon, prior.param, samples) {
   
   # first retrieve the final size of the observed data
   final.size.obs <- length(obs.data)
   
   # matrix to store the posterior samples
   post.samples <- matrix(NA, ncol = 2, nrow=samples)
   
   i = 0
   while (i < samples) {
      
      # draw from the prior distribution
      lambda <- rexp(1, prior.param[1])
      gamma <- runif(1, 0, 1)
      
      # simulate data
      sim.data <- simSIR.discrete(N, lambda, gamma, T)
      
      # get start/end dates
      start.date <- min(which(sim.data$pop[4,]==1)); start.date
      end.date <- start.date + len.out - 1; end.date
      
      # compute the distance
      d <- sum((obs.data - sim.data$pop[4,][start.date:end.date])^2)
      
      if (d < epsilon){
         i <- i + 1
         post.samples[i,] <- c(lambda, gamma)
            
      }
         
   }
   post.samples
}





##################################
# Coupled-ABC Homogeneously mixing SIR code
###################################

# Computes successive thresholds
thres=function(n)
{
   thres=rep(0,(n-1))
   thres[1]=rexp(1,(n-1)/n)
   for(i in 2:(n-1)) thres[i]=thres[i-1]+rexp(1,(n-i)/n)
   thres
}


 
# Simulates epidemics with a constant infectious period length 1.
# output gives the set of lambda parameters consistent with the data (m out of n infected)
# Only successful simulations kept 

epidemic=function(n,m,run)
{
   output=matrix(0,ncol=2,nrow=run)
   ss=seq(1,(n-1),1)
   count=0
   for(j in 1:run)
   {
      t=thres(n)
      y=t/ss
      q=max(y[1:(m-1)])
      if(q<y[m])
      {
         count=count+1
         output[count,]=c(q,y[m])
      }
   }
   print(count)
   output[1:count,]
}


 
# Simulates epidemics with a Gamma(k,k) infectious period.
# output gives the set of lambda parameters consistent with the data (m out of n infected)
# Only successful simulations kept 
epigamma=function(n,m,run,k)
{
   output=matrix(0,ncol=2,nrow=run)
   count=0
   for(j in 1:run)
   {
      t=thres(n)
      q=rgamma(n,k,k)
      y=0
      for(i in 1:(n-1)) y[i]=t[i]/sum(q[1:i])
      q=max(y[1:(m-1)])
      if(q<y[m])
      {
         count=count+1
         output[count,]=c(q,y[m])
      }
   }
   print(count)
   output[1:count,]
}


# Computing mean and variance

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




########### binning functions ######################

binning <- function (x, y, breaks, nbins) 
{
   binning.1d <- function(x, y, breaks, nbins) {
      f <- cut(x, breaks = breaks)
      if (any(is.na(f))) 
         stop("breaks do not span the range of x")
      freq <- tabulate(f, length(levels(f)))
      midpoints <- (breaks[-1] + breaks[-(nbins + 1)])/2
      id <- (freq > 0)
      x <- midpoints[id]
      x.freq <- as.vector(freq[id])
      result <- list(x = x, x.freq = x.freq, table.freq = freq, 
                     breaks = breaks)
      if (!all(is.na(y))) {
         result$means <- as.vector(tapply(y, f, mean))[id]
         result$sums <- as.vector(tapply(y, f, sum))[id]
         result$devs <- as.vector(tapply(y, f, function(x) sum((x - 
                                                                   mean(x))^2)))[id]
      }
      result
   }
   binning.2d <- function(x, y, breaks, nbins) {
      f1 <- cut(x[, 1], breaks = breaks[, 1])
      f2 <- cut(x[, 2], breaks = breaks[, 2])
      freq <- t(table(f1, f2))
      dimnames(freq) <- NULL
      midpoints <- (breaks[-1, ] + breaks[-(nbins + 1), ])/2
      z1 <- midpoints[, 1]
      z2 <- midpoints[, 2]
      X <- cbind(rep(z1, length(z2)), rep(z2, rep(length(z1), 
                                                  length(z2))))
      X.f <- as.vector(t(freq))
      id <- (X.f > 0)
      X <- X[id, ]
      dimnames(X) <- list(NULL, dimnames(x)[[2]])
      X.f <- X.f[id]
      result <- list(x = X, x.freq = X.f, midpoints = midpoints, 
                     breaks = breaks, table.freq = freq)
      if (!all(is.na(y))) {
         result$means <- as.numeric(tapply(y, list(f1, f2), 
                                           mean))[id]
         result$devs <- as.numeric(tapply(y, list(f1, f2), 
                                          function(x) sum((x - mean(x))^2)))[id]
      }
      result
   }
   binning.3d <- function(x, y, breaks, nbins) {
      f1 <- cut(x[, 1], breaks = breaks[, 1])
      f2 <- cut(x[, 2], breaks = breaks[, 2])
      f3 <- cut(x[, 3], breaks = breaks[, 3])
      freq <- table(f1, f2, f3)
      dimnames(freq) <- NULL
      midpoints <- (breaks[-1, ] + breaks[-(nbins + 1), ])/2
      z1 <- midpoints[, 1]
      z2 <- midpoints[, 2]
      z3 <- midpoints[, 3]
      X <- as.matrix(expand.grid(z1, z2, z3))
      X.f <- as.vector(freq)
      id <- (X.f > 0)
      X <- X[id, ]
      dimnames(X) <- list(NULL, dimnames(x)[[2]])
      X.f <- X.f[id]
      result <- list(x = X, x.freq = X.f, midpoints = midpoints, 
                     breaks = breaks, table.freq = freq)
      if (!all(is.na(y))) {
         result$means <- as.numeric(tapply(y, list(f1, f2, 
                                                   f3), mean))[id]
         result$devs <- as.numeric(tapply(y, list(f1, f2, 
                                                  f3), function(x) sum((x - mean(x))^2)))[id]
      }
      result
   }
   if (length(dim(x)) > 0) {
      if (!isMatrix(x)) 
         stop("wrong parameter x for binning")
      ndim <- dim(x)[2]
      if (ndim > 3) 
         stop("binning can be carried out only with 1-3 variables")
      if (missing(y)) 
         y <- rep(NA, nrow(x))
      if (missing(nbins)) 
         nbins <- round(log(nrow(x))/log(2) + 1)
      if (missing(breaks)) {
         breaks <- cbind(seq(min(x[, 1]), max(x[, 1]), length = nbins + 
                                1), seq(min(x[, 2]), max(x[, 2]), length = nbins + 
                                           1))
         if (ndim == 3) 
            breaks <- cbind(breaks, seq(min(x[, 3]), max(x[, 
                                                           3]), length = nbins + 1))
         breaks[1, ] <- breaks[1, ] - rep(10^(-5), ncol(breaks))
      }
      else nbins <- nrow(breaks) - 1
      if (max(abs(breaks)) == Inf | is.na(max(abs(breaks)))) 
         stop("illegal breaks")
      if (ndim == 2) 
         result <- binning.2d(x, y, breaks = breaks, nbins = nbins)
      else result <- binning.3d(x, y, breaks = breaks, nbins = nbins)
   }
   else {
      x <- as.vector(x)
      if (missing(y)) 
         y <- rep(NA, length(x))
      if (missing(nbins)) 
         nbins <- round(log(length(x))/log(2) + 1)
      if (missing(breaks)) {
         breaks <- seq(min(x), max(x), length = nbins + 1)
         breaks[1] <- breaks[1] - 10^(-5)
      }
      else nbins <- length(breaks) - 1
      if (max(abs(breaks)) == Inf | is.na(max(abs(breaks)))) 
         stop("illegal breaks")
      result <- binning.1d(x, y, breaks = breaks, nbins = nbins)
   }
   result
}
