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

##=============================================
# Initialize main variables and generate data #
##=============================================
K           <- 3      # Number of clusters
N           <- 500    # Number of objects
X   <- gen.gaussian(N=N, K=K, pi.c=c(.4,.3,.3), mus=c(0,2,5), stds=c(1,1,1))

##=========================
# Initialize parameters   #
##=========================
N.Sims          <- 10000        # Set the number of simulations
burnin          <- 5000         # Set how many samples should be discarded
dir.a           <- rep(1, K)    # Dirichlet concentration parameter
pi.cur          <- rep(1/K, K)  # Initialize mixing proportions for each cluster

Normal          <- list()       # Create an object of type Normal
Normal$mu       <- c(1,2,3)     # Initialize means for each cluster
for (k in 1:K){
  Normal$Tau[k] <- 1/sd(X)      # Normal precision for each cluster
}
Normal$Norm     <- list(mu.0=mean(X), tau.0=1/sd(X))  # Normal hyperparameters
Normal$Gamma    <- list(shape.0=1, rate.0=1)          # Gamma hyperparameters

params          <- list(Normal=Normal, pi.cur=pi.cur, dir.a=dir.a)

##=======================================================
# Do inference using Gibbs sampling, explicitly giving  #
# initial parameters.                                   #
##=======================================================
gibbs <- gmm1D.gibbs(X=X, 
                     K=K, 
                     N.Sims=N.Sims, 
                     burnin=burnin, 
                     params=params,
                     logl=TRUE)

##=================================================================
# Do inference using Gibbs sampling, without initial parameters,  #
# the method will initilize parameters using k-means              #
##=================================================================
gibbs.kmeans <- gmm1D.gibbs(X=X, K=K)


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
  density <- dnorm(x, mean=gibbs$summary$mu[k], sd=sqrt(1/gibbs$summary$tau[k]))
  mixture <- mixture + density * gibbs$summary$pi[k]
}
lines(x,mixture,col="red",lwd=2)


##===========================================
# Example of using cumsum function to plot  #
# the running mean of the MCMC samples.     #
##===========================================
plot(cumsum(gibbs$draws$mu[,1])/(1:length(gibbs$draws$mu[,1])), type="l", 
     xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(-1,10))
lines(cumsum(gibbs$draws$mu[,2])/(1:length(gibbs$draws$mu[,2])), lwd=2, col="orange3")
lines(cumsum(gibbs$draws$mu[,3])/(1:length(gibbs$draws$mu[,3])), lwd=2, col="darkgreen")

plot(gibbs$draws$pi[,1], type="l", xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(0,1))
lines(gibbs$draws$pi[,2], lwd=2, col="orange3")
lines(gibbs$draws$pi[,3], lwd=2, col="darkgreen")


##===========================================
# Create mcmc objects for the MCMC draws    #
# and use the coda for plotting and in      #
# general summary statistics                #
##===========================================
mu.draws <- mcmc(gibbs$draws$mu)
plot(mu.draws)
HPDinterval(mu.draws, 0.95)

tau.draws <- mcmc(gibbs$draws$tau)
plot(tau.draws)

pi.draws <- mcmc(gibbs$draws$pi)
plot(pi.draws)

combinedchains = mcmc.list(mcmc(gibbs$draws$mu), mcmc(gibbs.kmeans$draws$mu))
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)

NLL <- mcmc(gibbs$summary$NLL)
traceplot(NLL)

plotMixAutoc(gibbs)
