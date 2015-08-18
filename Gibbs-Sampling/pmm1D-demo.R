##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
library(coda)
library(label.switching)
library(R.utils)
source('pmm1D-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)
set.seed(12345)

##=============================================
# Initialize main variables and generate data #
##=============================================
K   <- 3        # Number of clusters
N   <- 1000     # Number of objects
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
ord.constr    <- TRUE
stephens      <- TRUE

##=======================================================
# Do inference using Gibbs sampling, explicitly giving  #
# initial parameters.                                   #
##=======================================================
gibbs <- pmm1D.gibbs(X=X, 
                     K=K, 
                     N.Sims=N.Sims, 
                     burnin=burnin, 
                     params=params,
                     ord.constr=ord.constr,
                     stephens=stephens,
                     logl=TRUE)

##=================================================================
# Do inference using Gibbs sampling, without initial parameters,  #
# the method will initilize parameters using k-means              #
##=================================================================
gibbs.kmeans <- pmm1D.gibbs(X=X, 
                            K=K, 
                            N.Sims=N.Sims, 
                            burnin=burnin, 
                            ord.constr=ord.constr,
                            stephens=stephens)

##=====================================
# Plot the data points and their pdfs #
##=====================================
invisible(readline(prompt="Press [enter] to show the plots"))
par()
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
     xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(0,20))
lines(cumsum(gibbs$draws$l[,2])/(1:length(gibbs$draws$l[,2])), lwd=2, col="orange3")
lines(cumsum(gibbs$draws$l[,3])/(1:length(gibbs$draws$l[,3])), lwd=2, col="darkgreen")

plot(gibbs$draws$pi[,1], type="l", xlab="time", ylab="x", lwd=2, col="steelblue", ylim=c(0,1))
lines(gibbs$draws$pi[,2], lwd=2, col="orange3")
lines(gibbs$draws$pi[,3], lwd=2, col="darkgreen")

##===================================================
# Plot autocorrelation functions using ACF method.  #
##===================================================
plotMixAutocorr(gibbs)

##===================================================
# Plot traces of sampled parameters during MCMC.    #
##===================================================
plotMixTrace(gibbs)

##=======================================================
# Create mcmc objects for the MCMC draws and use 'coda' #  
# for plotting and generate summary statistics          #
##=======================================================
l.draws <- mcmc(gibbs$draws$l)
plot(l.draws, main="Mean l")
HPDinterval(l.draws)

NLL <- mcmc(gibbs$summary$NLL)
traceplot(NLL, main="NLL trace plot")

##=================================================
# Use Gelman convergence diagnostic method. For   #
# this method we need more than one MCMC chains.  #      
##=================================================
if ((gibbs$dat$N.Sims - gibbs$dat$burnin) == (gibbs.kmeans$dat$N.Sims - gibbs.kmeans$dat$burnin)){
  combinedchains = mcmc.list(mcmc(gibbs$draws$l), mcmc(gibbs.kmeans$draws$l))
  plot(combinedchains)
  gelman.diag(combinedchains)
  gelman.plot(combinedchains)
}else{
  message("The MCMC chains have different length!")
}


##===========================================================
# Use 'label.switching' package to relabel the MCMC output  #
# due to identifiability issues when using Gibbs sampling   #
# for mixture models.                                       #
##===========================================================
mcmc.out <- array(0, dim=c((gibbs$dat$N.Sims-gibbs$dat$burnin), gibbs$dat$K, 2))
mcmc.out[,,1] <- gibbs$draws$pi
mcmc.out[,,2] <- gibbs$draws$l


## 1. Artificial Identifiability Constraints method
aic.relabel <- aic(mcmc = mcmc.out, constraint=2)
aic.mcmc    <- permute.mcmc(mcmc.out, aic.relabel$permutations)[[1]]

aic.pi  <- mcmc(aic.mcmc[,,1])
plot(aic.pi, main="AIC: Mix. prop. pi")
aic.l   <- mcmc(aic.mcmc[,,2])
plot(aic.l, main="AIC: Mean l")

if (stephens){
  ## 2. Stephens' algorithm using KL divergence
  ste.KL    <- stephens(gibbs$summary$PR)
  KL.mcmc   <- permute.mcmc(mcmc.out,ste.KL$permutations)[[1]]
  
  KL.pi <- mcmc(KL.mcmc[,,1])
  plot(KL.pi, main="KL: Mix. prop. pi")
  KL.l  <- mcmc(KL.mcmc[,,2])
  plot(KL.l, main="KL: Mean l")
}
