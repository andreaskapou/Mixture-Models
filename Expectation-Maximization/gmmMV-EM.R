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
pi      <- as.vector(table(C.n)/NROW(X)) # mixing proportions
mu      <- cl$centers                 # means for each Gaussian
Sigma   <- list(length=K)             # Covariance matrix for each Gaussian
for (k in 1:K){
  Sigma[[k]]  <- cov(X[C.n==k,])
}

##===============================
# Run Expectation Maximization  #
##===============================
# Hold responsibilities that component k takes on explaining the observation x.n
post.resp   <- matrix(, N, K)
log.likel   <- 0
ll          <- 0

for (i in 1:1000){  # Loop until convergence
  ##===============================================
  ## Expectation Step
  pdf.w     <- matrix(, N, K)   # Hold the PDF of each point on each cluster k
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- pi[k] * dmvnorm(X, mean=mu[k,], sigma=Sigma[[k]], log=F)
  }
  post.resp <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  
  ##===============================================
  ## Maximization Step
  prev.log.likel <- log.likel # Store to check for convergence
  for (k in 1:K){
    # Update mixing proportions
    pi[k]     <- mean(post.resp[, k])
    # Update mean for cluster 'k' by taking the weighted average of *all* data points.
    mu[k,]     <- (t(post.resp[,k]) %*% X) / sum(post.resp[,k])
    # Calculate the variance for cluster 'k'
    X.cen <- sweep(X, MARGIN=2, mu[k,], FUN="-")
    Sigma[[k]] <- t(X.cen) %*% (X.cen * post.resp[,k]) / sum(post.resp[,k])
    
    # Calculate the log likelihood
    ll[k] <- -.5 * sum( post.resp[,k] * dmvnorm(X, mu[k,], sigma=Sigma[[k]], log=T))
  }
  # Check for convergence.
  log.likel = sum(ll)
  if (abs(log.likel-prev.log.likel) < epsilon){
    break
  }
  print(log.likel)
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
    mixture <- mixture + pi[k]*dmvnorm(cbind(x,y),mean=mu[k,],sigma=Sigma[[k]])
  }
  return(mixture)
  })
contour(xpts,ypts,mixture.contour,nlevels=10,drawlabel=FALSE,col="red",
        xlab="Eruption time",ylab="Waiting time",
        main="Waiting vs Eruption time - Old Faithful")
points(faithful)