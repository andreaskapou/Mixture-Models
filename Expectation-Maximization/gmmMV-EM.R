##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
require(mvtnorm)
require(ggplot2)
library(MASS)
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory

##====================
# Generate the data  #
##====================
epsilon <- 0.00001    # Convergence paramater
X       <- as.matrix(faithful)
K       <- 2          # Number of clusters

##=========================
# Initialize variables    #
##=========================
N       <- nrow(X)                    # Length of the dataset
cl      <- kmeans(X, K, nstart = 25)  # Use Kmeans with random starts
C.n     <- cl$cluster                 # get the mixture components
pi.c    <- as.vector(table(C.n)/NROW(X)) # mixing proportions
mu      <- cl$centers                 # means for each Gaussian
Sigma   <- list(length=K)             # Covariance matrix for each Gaussian
for (k in 1:K){
  Sigma[[k]]  <- cov(X[C.n==k,])
}
post.resp   <- matrix(, N, K)         # Hold responsibilities
pdf.w       <- matrix(, N, K)         # Hold PDF of each point on each cluster k
logLik      <- 0                      # Initialize log likelihood

##===============================
# Run Expectation Maximization  #
##===============================
for (i in 1:1000){  # Loop until convergence
  prevLogLik  <- logLik               # Store to check for convergence
  
  ##========
  # E-Step #
  ##========
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- pi.c[k] * dmvnorm(X, mean=mu[k,], sigma=Sigma[[k]], log=F)
  }
  post.resp <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  
  ##========
  # M-Step #
  ##========
  for (k in 1:K){
    N.k       <- sum(post.resp[,k])     # Sum of responsibilities for cluster k
    pi.c[k]   <- N.k / N                # Update mixing proportions for cluster k
    mu[k,]     <- (t(post.resp[,k]) %*% X) / sum(post.resp[,k]) # Update mu
    X.cen <- sweep(X, MARGIN=2, mu[k,], FUN="-") # Update Sigma
    Sigma[[k]] <- t(X.cen) %*% (X.cen * post.resp[,k]) / sum(post.resp[,k])
  }
  
  # Evaluate the log likelihood
  # ln p(X|mu,S,p) = Sum_{n=1}^{N}(ln(Sum_{k=1}^{K}(p_k * N(x_n|mu_k, S_k))))
  logLik   <- sum(log(colSums(pdf.w)))
  if (abs(logLik-prevLogLik) < epsilon){ # Check for convergence.
    break
  }
  print(logLik)
} #End of Expectation Maximization loop.

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
    mixture <- mixture + pi.c[k]*dmvnorm(cbind(x,y),mean=mu[k,],sigma=Sigma[[k]])
  }
  return(mixture)
  })
contour(xpts,ypts,mixture.contour,nlevels=10,drawlabel=FALSE,col="red",
        xlab="Eruption time",ylab="Waiting time",
        main="Waiting vs Eruption time - Old Faithful")
points(faithful)