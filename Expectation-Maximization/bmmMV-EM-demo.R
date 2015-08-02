##=========================================================
# Script for running Binomial Mixture Models using EM to  #
# fit the model to the data.                              #
##=========================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("bmmMV-EM.R")
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(123456)         # Set seed for reproducible results

##=============================================
# Generate the data and initialize variables #
##=============================================
N       <- 1000         # Number of data points
K       <- 2            # Number of clusters
r1      <- rbinom(n=N, size=40, prob=.9)  # Total number of trials 1st dim
r2      <- rbinom(n=N, size=60, prob=.9)  # Total number of trials 2nd dim
r       <- as.matrix(cbind(r1, r2))
X       <- gen.MV.binomial(N=N, r=r)

epsilon <- 1e-10        # Convergence paramater for EM
maxIter <- 1000         # Maximum number of iterations for EM

# Optional parameter initialization.
p       <- matrix(c(.1,.2,.3,.4), ncol=2)  # Mean of each cluster
pi.c    <- rep(1/K, K)  # Initialize mixing proportions for each cluster

theta   <- list(p=p, pi.c=pi.c) # Wrap all the parameters in a list

##===================================================
# Run BMM-EM, explicitly giving initial parameters  #
##===================================================
fit.bmm <- bmmMV.EM(X=X,
                  r=r,
                  K=K, 
                  theta=theta, 
                  epsilon=epsilon, 
                  maxIter=maxIter, 
                  isDebug=TRUE)

##===================================================
# Run GMM-EM, without initial parameters, then the  #
# method will initilize parameters using k-means    #
##===================================================
fit.bmm.kmeans <- bmmMV.EM(X=X, r=r, K=K)


cat("\n\nBMM-EM (with given parameters):\n")
cat("Total iterations:", length(fit.bmm$all.NLL), "\n")
cat("Probabilities:\n")
cat("     Cluster 1:", fit.bmm$p[1,], "\n")
cat("     Cluster 2:", fit.bmm$p[2,], "\n")
cat("Mixing proportions:", fit.bmm$pi.c, "\n")
cat("Best NLL:", fit.bmm$NLL, "\n")
cat("BIC:", fit.bmm$BIC, "\n\n")

cat("BMM-EM (with parameters initialized by k-means):\n")
cat("Total iterations:", length(fit.bmm.kmeans$all.NLL), "\n")
cat("Probabilities:\n")
cat("     Cluster 1:", fit.bmm$p[1,], "\n")
cat("     Cluster 2:", fit.bmm$p[2,], "\n")
cat("Mixing proportions:", fit.bmm.kmeans$pi.c, "\n")
cat("Best NLL:", fit.bmm.kmeans$NLL, "\n")
cat("BIC:", fit.bmm.kmeans$BIC, "\n\n")
