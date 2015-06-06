##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
library(coda)
library(R.utils)
source('gmm1D-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)

set.seed(1)
##===========================
# Initialize main variables #
##===========================
K           <- 3      # Number of clusters
N           <- 500    # Number of objects
N.Sims      <- 10000  # Set the number of simulations
burnin      <- 5000   # Set how many samples should be burned in
Normal      <- list() # Create a Normal object

##====================
# Generate the data  #
##====================
X <- gen.gaussian(N=N, K=K, pi.c=c(.4,.3,.3), mus=c(0,2,5), stds=c(1,1,1))

##=========================
# Initialize parameters   #
##=========================
cl              <- kmeans(X, K, nstart = 25)    # Use Kmeans with random starts
C.n             <- cl$cluster                   # get the mixture components
pi.cur          <- as.vector(table(C.n)/NROW(X))# mixing proportions
dir.a           <- rep(1/K, K)                  # Dirichlet concentration parameter
Normal$mu       <- as.vector(cl$centers)        # Normal mean for each cluster
for (k in 1:K){
  Normal$Tau[k] <- 1/var(X[C.n==k])             # Normal precision for each cluster
}
Normal$Norm     <- list(mu.0=mean(X), tau.0=1/sd(X))  # Normal hyperparameters
Normal$Gamma    <- list(a=1, b=1)               # Gamma hyperparameters
logl            <- TRUE                         # If we want to compute log likel

##===================================
# Do inference using Gibbs sampling #
##===================================
gibbs <- gmm1D.gibbs(X, K, N.Sims, burnin, Normal, pi.cur, dir.a, logl=logl)

##=====================================
# Plot the data points and their pdfs #
##=====================================
invisible(readline(prompt="Press [enter] to show the plots"))
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 0.1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0.0, .20), main=paste("Density plot of",K,"Gaussian Mixtures"),xlab="x")
mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dnorm(x, mean=mean(gibbs$mu.draws[,k]),sd=sqrt(1/mean(gibbs$tau.draws[,k])))
  mixture <- mixture + density * mean(gibbs$pi.draws[,k])
}
lines(x,mixture,col="red",lwd=2)


##===========================================
# Example of using cumsum function to plot  #
# the running mean of the MCMC samples.     #
##===========================================
plot(cumsum(gibbs$mu.draws[,1])/(1:length(gibbs$mu.draws[,1])), type="l", 
     xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(-1,9))
lines(cumsum(gibbs$mu.draws[,2])/(1:length(gibbs$mu.draws[,2])), lwd=2, col="orange3")
lines(cumsum(gibbs$mu.draws[,3])/(1:length(gibbs$mu.draws[,3])), lwd=2, col="darkgreen")

plot(gibbs$pi.draws[,1], type="l", xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(0.2,0.5))
lines(gibbs$pi.draws[,2], lwd=2, col="orange3")
lines(gibbs$pi.draws[,3], lwd=2, col="darkgreen")


##===========================================
# Create mcmc objects for the MCMC draws    #
# and use the coda for plotting and in      #
# general summary statistics                #
##===========================================
mu.draws <- mcmc(gibbs$mu.draws)
plot(mu.draws)

tau.draws <- mcmc(gibbs$tau.draws)
plot(tau.draws)

pi.draws <- mcmc(gibbs$pi.draws)
plot(pi.draws)

#chainmu4 <- mu.draws
#chaintau4 <- tau.draws
#chainpi4 <- pi.draws

#combinedchains = mcmc.list(chainmu1, chainmu2, chainmu3, chainmu4)
#plot(combinedchains)
#gelman.diag(combinedchains)
#gelman.plot(combinedchains)
