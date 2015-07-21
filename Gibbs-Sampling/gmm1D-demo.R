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
source('gmm1D-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)
set.seed(12345)

##=============================================
# Initialize main variables and generate data #
##=============================================
K   <- 3      # Number of clusters
N   <- 1000   # Number of objects
X   <- gen.gaussian(N=N, K=K, pi.c=c(.4,.3,.3), mus=c(0,5,10), stds=c(1,2,2))

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
ord.constr      <- FALSE
stephens        <- TRUE

##=======================================================
# Do inference using Gibbs sampling, explicitly giving  #
# initial parameters.                                   #
##=======================================================
gibbs <- gmm1D.gibbs(X=X, 
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
gibbs.kmeans <- gmm1D.gibbs(X=X, 
                            K=K, 
                            N.Sims=N.Sims, 
                            burnin=burnin, 
                            ord.constr=ord.constr,
                            stephens=stephens)



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
mu.draws <- mcmc(gibbs$draws$mu)
traceplot(mu.draws, main="Mean mu")
HPDinterval(mu.draws)

NLL <- mcmc(gibbs$summary$NLL)
traceplot(NLL, main="NLL trace plot")

##=================================================
# Use Gelman convergence diagnostic method. For   #
# this method we need more than one MCMC chains.  #      
##=================================================
if ((gibbs$dat$N.Sims - gibbs$dat$burnin) == (gibbs.kmeans$dat$N.Sims - gibbs.kmeans$dat$burnin)){
  # Relabel the cluster assignments so that we have label agreement between the different chains
  chain2 <- gibbs.kmeans$draws$mu
  chain2[,1] <- gibbs.kmeans$draws$mu[,3]
  chain2[,3] <- gibbs.kmeans$draws$mu[,1]
  combinedchains = mcmc.list(mcmc(gibbs$draws$mu), mcmc(chain2))
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
mcmc.out <- array(0, dim=c((gibbs$dat$N.Sims-gibbs$dat$burnin), gibbs$dat$K, 3))
mcmc.out[,,1] <- gibbs$draws$pi
mcmc.out[,,2] <- gibbs$draws$mu
mcmc.out[,,3] <- gibbs$draws$tau


## 1. Artificial Identifiability Constraints method
aic.relabel <- aic(mcmc = mcmc.out, constraint=2)
aic.mcmc    <- permute.mcmc(mcmc.out, aic.relabel$permutations)[[1]]

aic.pi  <- mcmc(aic.mcmc[,,1])
traceplot(aic.pi, main="AIC: Mix. prop. pi")
aic.mu  <- mcmc(aic.mcmc[,,2])
traceplot(aic.mu, main="AIC: Mean mu")
aic.tau <- mcmc(aic.mcmc[,,3])
traceplot(aic.tau, main="AIC: Precision tau")

if (stephens){
  ## 2. Stephens' algorithm using KL divergence
  ste.KL    <- stephens(gibbs$summary$PR)
  KL.mcmc   <- permute.mcmc(mcmc.out,ste.KL$permutations)[[1]]
  
  KL.pi <- mcmc(KL.mcmc[,,1])
  traceplot(KL.pi, main="KL: Mix. prop. pi")
  KL.mu <- mcmc(KL.mcmc[,,2])
  traceplot(KL.mu, main="KL: Mean mu")
  KL.tau <- mcmc(KL.mcmc[,,3])
  traceplot(KL.tau, main="KL: Precision tau")
}
