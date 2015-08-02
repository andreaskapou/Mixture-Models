##=========================================================
# Script for running Poisson Mixture Models using EM to   #
# fit the model to the data.                              #
##=========================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("pmmMV-EM.R")
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)         # Set seed for reproducible results

##=============================================
# Generate the data and initialize variables  #
##=============================================
N       <- 1000         # Number of data points
K       <- 2            # Number of clusters
X       <- gen.MV.poisson(N=N)

epsilon <- 1e-10        # Convergence paramater for EM
maxIter <- 1000         # Maximum number of iterations for EM

# Optional parameter initialization.
lambdas <- matrix(c(1,2,3,4),ncol=2)     # Mean and variance of each cluster
pi.c    <- rep(1/K, K)  # Initialize mixing proportions for each cluster

theta   <- list(lambdas=lambdas, pi.c=pi.c) # Wrap all the parameters in a list

##===================================================
# Run PMM-EM, explicitly giving initial parameters  #
##===================================================
fit.pmm <- pmmMV.EM(X=X,
                    K=K, 
                    theta=theta, 
                    epsilon=epsilon, 
                    maxIter=maxIter,  
                    isDebug=TRUE)

##===================================================
# Run PMM-EM, without initial parameters, then the  #
# method will initilize parameters using k-means    #
##===================================================
fit.pmm.kmeans <- pmmMV.EM(X=X, K=K)


cat("\n\nPMM-EM (with given parameters):\n")
cat("Total iterations:", length(fit.pmm$all.NLL), "\n")
cat("Mean:\n")
cat("     Cluster 1:", fit.pmm$lambdas[1,], "\n")
cat("     Cluster 2:", fit.pmm$lambdas[2,], "\n")
cat("Mixing proportions:", fit.pmm$pi.c, "\n")
cat("Best NLL:", fit.pmm$NLL, "\n")
cat("BIC:", fit.pmm$BIC, "\n\n")


cat("\n\nPMM-EM (with given parameters):\n")
cat("Total iterations:", length(fit.pmm.kmeans$all.NLL), "\n")
cat("Mean:\n")
cat("     Cluster 1:", fit.pmm.kmeans$lambdas[1,], "\n")
cat("     Cluster 2:", fit.pmm.kmeans$lambdas[2,], "\n")
cat("Mixing proportions:", fit.pmm.kmeans$pi.c, "\n")
cat("Best NLL:", fit.pmm.kmeans$NLL, "\n")
cat("BIC:", fit.pmm.kmeans$BIC, "\n\n")
