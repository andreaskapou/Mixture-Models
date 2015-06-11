##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
library(coda)
library(R.utils)
source('pmm1D-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)
set.seed(1235)

##=============================================
# Initialize main variables and generate data #
##=============================================
K           <- 3      # Number of clusters
N           <- 500    # Number of objects
X   <- gen.poisson(N=N, K=K, pi.c=c(.3,.2,.5), lambdas=c(4,10,15))

##=========================
# Initialize variables    #
##=========================
N.Sims        <- 10000        # Set the number of simulations
burnin        <- 5000         # Set how many samples should be discarded
dir.a         <- rep(1, K)    # Dirichlet concentration parameter
pi.cur        <- rep(1/K, K)  # Initialize mixing proportions for each cluster


Poisson       <- list()       # Create an object of type Poisson
Poisson$l     <- c(2,4,6)     # Initialize means for each cluster
Poisson$Gamma <- list(shape.0=1, rate.0=0.01)  # Gamma hyperparameters

params        <- list(Poisson=Poisson, pi.cur=pi.cur, dir.a=dir.a)

##=======================================================
# Do inference using Gibbs sampling, explicitly giving  #
# initial parameters.                                   #
##=======================================================
gibbs <- pmm1D.gibbs(X=X, 
                     K=K, 
                     N.Sims=N.Sims, 
                     burnin=burnin, 
                     params=params,
                     logl=TRUE)

##=================================================================
# Do inference using Gibbs sampling, without initial parameters,  #
# the method will initilize parameters using k-means              #
##=================================================================
gibbs.kmeans <- pmm1D.gibbs(X=X, K=K)

##=====================================
# Plot the data points and their pdfs #
##=====================================
invisible(readline(prompt="Press [enter] to show the plot"))
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0.0, .10), main=paste("Density plot of",K,"Poisson Mixtures"),xlab="x")

mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dpois(x, lambda=gibbs$summary$l[k])
  mixture <- mixture + density * gibbs$summary$pi[k]
}
lines(x,mixture,col="red",lwd=2)


##===========================================
# Example of using cumsum function to plot  #
# the running mean of the MCMC samples.     #
##===========================================
plot(cumsum(gibbs$draws$l[,1])/(1:length(gibbs$draws$l[,1])), type="l", 
     xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(-1,9))
lines(cumsum(gibbs$draws$l[,2])/(1:length(gibbs$draws$l[,2])), lwd=2, col="orange3")
lines(cumsum(gibbs$draws$l[,3])/(1:length(gibbs$draws$l[,3])), lwd=2, col="darkgreen")

plot(gibbs$draws$pi[,1], type="l", xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(0.1,0.6))
lines(gibbs$draws$pi[,2], lwd=2, col="orange3")
lines(gibbs$draws$pi[,3], lwd=2, col="darkgreen")


##===========================================
# Create mcmc objects for the MCMC draws    #
# and use the coda for plotting and in      #
# general summary statistics                #
##===========================================
lambda.draws <- mcmc(gibbs$draws$l)
plot(lambda.draws)
HPDinterval(lambda.draws, 0.95)

pi.draws <- mcmc(gibbs$draws$pi)
plot(pi.draws)

combinedchains = mcmc.list(mcmc(gibbs$draws$l), mcmc(gibbs.kmeans$draws$l))
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)

NLL <- mcmc(gibbs$summary$NLL)
traceplot(NLL)

plotMixAutoc(gibbs)
