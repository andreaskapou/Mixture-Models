##===================================================================
# Script for running Gaussian Mixture Models using EM to fit the    #
# model to the data. This script is only for multivariate Gaussian  #
# distributions, check the gmm-EM-demo.R script for univariate      #
# Gaussian distributions. At the end of the script, the results     #
# are also compared with the 'mixtools' package.                    #
##===================================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("gmmMV-EM.R")
library(R.utils)
library(mixtools)
library(mvtnorm)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)         # Set seed for reproducible results

X       <- as.matrix(faithful)
K       <- 2          # Number of clusters

##=====================================================
# Run GMM-MV-EM, without initial parameters, then the #
# method will initilize parameters using k-means      #
##=====================================================
fit.gmmMV.kmeans <- gmmMV.EM(X=X, K=K, isLog=FALSE, isDebug=TRUE)

##=====================================
# Plot the data points and their pdfs #
##=====================================
# Points to estimate the density and plot the contours
xpts <- seq(from=min(X[,1])-1,to=max(X[,1])+1,length.out=100)
ypts <- seq(from=min(X[,2])-1,to=max(X[,2])+2,length.out=100)

# Estimate the mixture contours
mixture.contour <- outer(xpts,ypts,function(x,y) {
  mixture <- matrix(data=0, length(xpts)*length(ypts), 1)
  for(k in 1:K){
    mixture <- mixture + fit.gmmMV.kmeans$pi.c[k]*mvtnorm::dmvnorm(cbind(x,y),
                                         mean=fit.gmmMV.kmeans$mu[k,],
                                         sigma=fit.gmmMV.kmeans$Sigma[[k]])
  }
  return(mixture)
})
contour(xpts,ypts,mixture.contour,nlevels=10,drawlabel=FALSE,col="red",
        xlab="Eruption time",ylab="Waiting time",
        main="Waiting vs Eruption time - Old Faithful")
points(faithful)


##===========================================
# Compare results with the mixtools package #
##===========================================
gm <- mvnormalmixEM(X, k=K)


cat("My GMM-EM (with parameters initialized by k-means):\n")
cat("Total iterations:", length(fit.gmmMV.kmeans$all.NLL), "\n")
cat("Mean:", fit.gmmMV.kmeans$mu, "\n")
#cat("Covariance:", fit.gmmMV.kmeans$s2, "\n")
cat("Mixing proportions:", fit.gmmMV.kmeans$pi.c, "\n")
cat("Best NLL:", fit.gmmMV.kmeans$NLL, "\n\n")

cat("Mixtools package GMM-EM:\n")
cat("Total iterations:", length(gm$all.loglik),"\n")
cat("Mean:", gm$mu[[1]], gm$mu[[1]], "\n")
#cat("Coavariance:", gm$sigma, "\n")
cat("Mixing proportions:", gm$lambda, "\n")
cat("Best NLL:", -gm$loglik, "\n\n")