##=========================================================
# Script for running Binomial Mixture Models using EM to  #
# fit the model to the data.                              #
##=========================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("bmm-EM.R")
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)         # Set seed for reproducible results

##=============================================
# Generate the data and initialize variables #
##=============================================
N       <- 500          # Number of data points
K       <- 3            # Number of clusters
r1      <- rbinom(n=N, size=40, prob=.9)  # Total number of trials
r       <- matrix(r1, ncol=1)
X       <- gen.binomial(N=N, K=3, pi.c=c(.4,.3,.3), p=c(0.2,0.7, 0.5), r=r)

epsilon <- 1e-10        # Convergence paramater for EM
maxIter <- 1000         # Maximum number of iterations for EM

# Optional parameter initialization.
p       <- c(.1,.2,.3)  # Mean of each cluster
pi.c    <- rep(1/K, K)  # Initialize mixing proportions for each cluster

theta   <- list(p=p, pi.c=pi.c) # Wrap all the parameters in a list

##===================================================
# Run BMM-EM, explicitly giving initial parameters  #
##===================================================
fit.bmm <- bmm.EM(X=X,
                  r=r,
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
fit.bmm.kmeans <- bmm.EM(X=X, r=r, K=K)


cat("\n\nBMM-EM (with given parameters):\n")
cat("Total iterations:", length(fit.bmm$all.NLL), "\n")
cat("Probabilities:", fit.bmm$mu, "\n")
cat("Mixing proportions:", fit.bmm$pi.c, "\n")
cat("Best NLL:", fit.bmm$NLL, "\n\n")

cat("BMM-EM (with parameters initialized by k-means):\n")
cat("Total iterations:", length(fit.bmm.kmeans$all.NLL), "\n")
cat("Probabilities:", fit.bmm.kmeans$pi, "\n")
cat("Mixing proportions:", fit.bmm.kmeans$pi.c, "\n")
cat("Best NLL:", fit.bmm.kmeans$NLL, "\n\n")
