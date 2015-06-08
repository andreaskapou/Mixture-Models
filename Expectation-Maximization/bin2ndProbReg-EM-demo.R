##=====================================================================
# Script for running Binomial Distributed Probit Regression Mixture   #
# Models using EM to fit the model to the data. This script runs only #
# for 2nd order polynomial functions. The general case for nth order  #
# polynomial functions can be run using binProbReg-EM-demo.R script   #
##=====================================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("bin2ndProbReg-EM.R")
library(R.utils)
library(ggplot2)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)                         # Set seed for reproducible results

##=============================================
# Generate the data  and initialize variables #
##=============================================
N       <- 400                          # Number of data points
K       <- 3                            # Number of clusters
X       <- gen.meth.data(N=N, pi.c=c(.5,.3,.2))  # Generate methylation profiles

epsilon <- 1e-4                         # Convergence paramater for EM
maxIter <- 1000                         # Maximum number of iterations for EM

# Optional parameter initialization.
pi.c    <- rep(1/K, K)                  # Initialize mixing proportions for each cluster
theta   <- matrix(0, nrow=3, ncol=K)    # Parameters of the 2nd order polynomial
for (k in 1:K){
  theta[,k] <- c(0+k/100, 0+k/100, 0+k/100)
}

params <- list(pi.c=pi.c, theta=theta)  # Wrap all the parameters in a list

##===========================================================
# Run BIN.PROB.REG-EM, explicitly giving initial parameters #
##===========================================================
fit.bin2ndProbReg <- bin2ndProbReg.EM(X=X, 
                                      K=K, 
                                      params=params, 
                                      epsilon=epsilon, 
                                      maxIter=maxIter,
                                      isDebug=TRUE)


##======================================================================
# Plot the results showing the K different functions that were learned #
##======================================================================
xs <- seq(-1,1,len=2000) # create some values
# Add extra space to right of plot area; change clipping to figure
#par(xpd=T, mar=par()$mar+c(0,0,0,4))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(X[[7]][1,], X[[7]][3,]/X[[7]][2,], col="darkgreen", pch=24, xlim=c(-1,1), 
     ylim=c(0.05,0.92), xlab="region x", ylab="methylation level")
lines(x=xs, y=probPolynomFun(theta=fit.bin2ndProbReg$theta[,1], xs), 
      col="darkgreen", lwd=2)
points(X[[245]][1,], X[[245]][3,]/X[[245]][2,], pch=23, col="red")
lines(x=xs, y=probPolynomFun(theta=fit.bin2ndProbReg$theta[,3], xs), 
      col="red", lwd=2)
points(X[[390]][1,], X[[390]][3,]/X[[390]][2,], pch=8, col="darkblue")
lines(x=xs, y=probPolynomFun(theta=fit.bin2ndProbReg$theta[,2], xs), 
      col="darkblue", lwd=2)

# Add legend to top right, outside plot region
#legend(locator(1),c("Cluster", "group B"), pch = c(1,2), lty = c(1,2))
legend("right", inset=c(-0.1,0), legend=c("1","2", "3"), pch=c(24,23,8), 
       col=c("darkgreen", "red", "darkblue"), title="Cluster")


##=============================================================
# Store plots and results locally in the disk for future use  #
##=============================================================
dev.copy(png, width = 800, height = 600, 
         filename=paste("../images/probit2ndParams", 
                        format(Sys.time(), "%a%b%d%H%M"),".png", sep=""))
dev.off()
save(fit.bin2ndProbReg, file=paste("../files/probit2ndTheta", 
                       format(Sys.time(), "%a%b%d%H%M"),".RData", sep=""))
