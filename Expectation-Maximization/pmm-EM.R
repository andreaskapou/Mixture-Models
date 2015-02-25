##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("../readData.R")

##====================
# Generate the data  #
##====================
epsilon <- 0.00001    # Convergence paramater
K       <- 3          # Number of clusters
X       <- gen.poisson(K=K, pi.c=c(.3,.2,.5), lambdas=c(5,20,50))

##=========================
# Initialize variables    #
##=========================
N       <- length(X)  # Length of the dataset
cl      <- kmeans(X, K, nstart = 25) # Use Kmeans with random starts
C.n     <- cl$cluster   # get the mixture components
pi.c    <- as.vector(table(C.n)/length(X)) # mixing proportions
lambdas <- as.vector(cl$centers) # means for each Gaussian

##===============================
# Run Expectation Maximization  #
##===============================
# Hold responsibilities that component k takes on explaining the observation x.n
post.resp   <- matrix(, N, K)
log.likel   <- 0

for (i in 1:1000){  # Loop until convergence
  ##===============================================
  ## Expectation step
  pdf.w     <- matrix(, N, K)   # Hold the PDF of each point on each cluster k
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- pi.c[k] * dpois(X, lambda=lambdas[k])
  }
  post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  
  ##===============================================
  ## Maximization step
  prev.log.likel <- log.likel # Store to check for convergence
  for (k in 1:K){
    # Update mixing proportions
    pi.c[k]     <- sum(post.resp[, k] / sum(post.resp))
    # Update mean for cluster 'k' by taking the weighted average of *all* data points.
    val         <- t(post.resp[, k]) %*% X
    lambdas[k]  <- val / sum(post.resp[,k])
  }
  # Check for convergence.
  log.likel   <- sum(log(sum(pdf.w)))
  if (abs(log.likel-prev.log.likel) < epsilon){
    break
  }
  print(log.likel)
} #End of Expectation Maximization loop.

##=====================================
# Plot the data points and their pdfs #
##=====================================
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0, .15), main="Density plot of 2 Gaussians", xlab = "x")

mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dpois(x, lambda=lambdas[k])
  mixture <- mixture + density * pi.c[k]
}
lines(x, mixture, col="red",lwd=2)