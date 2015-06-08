##=========================================================
# Script for running Poisson Mixture Models using EM to   #
# fit the model to the data.                              #
##=========================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("pmm-EM.R")
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)         # Set seed for reproducible results

##=============================================
# Generate the data and initialize variables  #
##=============================================
N       <- 500          # Number of data points
K       <- 3            # Number of clusters
X       <- gen.poisson(K=K, pi.c=c(.3,.2,.5), lambdas=c(7,15,28))

epsilon <- 1e-10        # Convergence paramater for EM
maxIter <- 1000         # Maximum number of iterations for EM

# Optional parameter initialization.
lambdas <- c(1,2,3)     # Mean and variance of each cluster
pi.c    <- rep(1/K, K)  # Initialize mixing proportions for each cluster

theta   <- list(lambdas=lambdas, pi.c=pi.c) # Wrap all the parameters in a list

##===================================================
# Run PMM-EM, explicitly giving initial parameters  #
##===================================================
fit.pmm <- pmm.EM(X=X,
                  K=K, 
                  theta=theta, 
                  epsilon=epsilon, 
                  maxIter=maxIter, 
                  isLog=TRUE, 
                  isDebug=TRUE)

##===================================================
# Run PMM-EM, without initial parameters, then the  #
# method will initilize parameters using k-means    #
##===================================================
fit.pmm.kmeans <- pmm.EM(X=X, K=K)


##=====================================
# Plot the data points and their pdfs #
##=====================================
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0, .15), main=paste("Density plot of",K,"Poissons"), xlab = "x")

mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dpois(x, lambda=fit.pmm$lambdas[k])
  mixture <- mixture + density * fit.pmm$pi.c[k]
}
lines(x, mixture, col="red",lwd=2)


cat("\n\nPMM-EM (with given parameters):\n")
cat("Total iterations:", length(fit.pmm$all.NLL), "\n")
cat("Mean:", fit.pmm$lambdas, "\n")
cat("Mixing proportions:", fit.pmm$pi.c, "\n")
cat("Best NLL:", fit.pmm$NLL, "\n\n")

cat("PMM-EM (with parameters initialized by k-means):\n")
cat("Total iterations:", length(fit.pmm.kmeans$all.NLL), "\n")
cat("Mean:", fit.pmm.kmeans$lambdas, "\n")
cat("Mixing proportions:", fit.pmm.kmeans$pi.c, "\n")
cat("Best NLL:", fit.pmm.kmeans$NLL, "\n\n")
