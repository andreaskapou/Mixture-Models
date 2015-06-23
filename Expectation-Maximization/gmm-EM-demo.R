##===================================================================
# Script for running Gaussian Mixture Models using EM to fit the    #
# model to the data. This script is only for univariate Gaussian    #
# distributions, check the gmmMV-EM-demo.R script for multivariate  #
# Gaussian distributions. At the end of the script, the results     #
# are also compared with the 'mixtools' package.                    #
##===================================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("gmm-EM.R")
library(R.utils)
library(mixtools)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)         # Set seed for reproducible results

##=============================================
# Generate the data and initialize variables #
##=============================================
N       <- 500          # Number of data points
K       <- 3            # Number of clusters
X       <- gen.gaussian(N=N, K=K, pi.c=c(.4,.3,.3), mus=c(-1,3,6), stds=c(1,1.1,1.2))

epsilon <- 1e-10        # Convergence paramater for EM
maxIter <- 1000         # Maximum number of iterations for EM

# Optional parameter initialization.
mu      <- c(1,2,3)     # Initialize means for each cluster
s2      <- c(2,2,2)     # Initialize variances for each cluster
pi.c    <- rep(1/K, K)  # Initialize mixing proportions for each cluster

theta   <- list(mu=mu, s2=s2, pi.c=pi.c) # Wrap all the parameters in a list

##===================================================
# Run GMM-EM, explicitly giving initial parameters  #
##===================================================
fit.gmm <- gmm.EM(X=X, 
                  K=K, 
                  theta=theta, 
                  epsilon=epsilon, 
                  maxIter=maxIter, 
                  isLog=TRUE, 
                  isDebug=TRUE)

##===================================================
# Run GMM-EM, without initial parameters, then the  #
# method will initilize parameters using k-means    #
##===================================================
fit.gmm.kmeans <- gmm.EM(X=X, K=K)


##=====================================
# Plot the data points and their pdfs #
##=====================================
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 0.1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0.0, .20), main=paste("Density plot of",K,"GMMs"), xlab = "x")
mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dnorm(x, mean=fit.gmm$mu[k],sd=sqrt(fit.gmm$s2[k]))
  mixture <- mixture + density * fit.gmm$pi.c[k]
}
lines(x,mixture,col="red",lwd=2)


##===========================================
# Compare results with the mixtools package #
##===========================================
gm <- normalmixEM(X, 
                  k=K, 
                  lambda=pi.c, 
                  mu=mu, 
                  sigma=sqrt(s2))

cat("\n\nMy GMM-EM (with given parameters):\n")
cat("Total iterations:", length(fit.gmm$all.NLL), "\n")
cat("Mean:", fit.gmm$mu, "\n")
cat("Variance:", fit.gmm$s2, "\n")
cat("Mixing proportions:", fit.gmm$pi.c, "\n")
cat("Best NLL:", fit.gmm$NLL, "\n")
cat("BIC:", fit.gmm$BIC, "\n\n")

cat("My GMM-EM (with parameters initialized by k-means):\n")
cat("Total iterations:", length(fit.gmm.kmeans$all.NLL), "\n")
cat("Mean:", fit.gmm.kmeans$mu, "\n")
cat("Variance:", fit.gmm.kmeans$s2, "\n")
cat("Mixing proportions:", fit.gmm.kmeans$pi.c, "\n")
cat("Best NLL:", fit.gmm.kmeans$NLL, "\n")
cat("BIC:", fit.gmm.kmeans$BIC, "\n\n")

cat("Mixtools package GMM-EM:\n")
cat("Total iterations:", length(gm$all.loglik),"\n")
cat("Mean:", gm$mu, "\n")
cat("Variance:", gm$sigma, "\n")
cat("Mixing proportions:", gm$lambda, "\n")
cat("Best NLL:", -gm$loglik, "\n\n")
