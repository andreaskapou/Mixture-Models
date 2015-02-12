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
X       <- gen.gaussian(K=K, pi=c(.4,.3,.3), mus=c(10,20,30), stds=c(2,3,1))

##=========================
# Initialize variables    #
##=========================
N       <- length(X)                  # Length of the dataset
cl      <- kmeans(X, K, nstart = 25)  # Use Kmeans with random starts
C.n     <- cl$cluster                 # get the mixture components
mu      <- as.vector(cl$centers)      # means for each Gaussian
pi      <- as.vector(table(C.n)/length(X)) # mixing proportions
sigma2  <- vector(length=K)           # Variance for each Gaussian
for (k in 1:K){
  sigma2[k]  <- var(X[C.n==k])
}
##===============================
# Run Expectation Maximization  #
##===============================
# Hold responsibilities that component k takes on explaining the observation x.n
post.resp   <- matrix(, N, K)
log.likel   <- 0

for (i in 1:1000){  # Loop until convergence
  ##===============================================
  ## Expectation Step
  pdf.w       <- matrix(, N, K)   # Hold the PDF of each point on each cluster k
  for (k in 1:K){ # Calculate the weighted PDF of each cluster for each data point
    pdf.w[,k] <- pi[k] * dnorm(X, mean=mu[k], sd=sqrt(sigma2[k]))
  }
  post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  
  ##===============================================
  ## Maximization Step
  prev.log.likel <- log.likel # Store to check for convergence
  for (k in 1:K){
    # Update mixing proportions
    pi[k]     <- mean(post.resp[, k])
    # Update mean for cluster 'k' by taking the weighted average of *all* data points.
    val       <- t(post.resp[, k]) %*% X
    mu[k]     <- val / sum(post.resp[, k])
    # Calculate the variance for cluster 'k' by taking the weighted
    # average of the squared differences from the mean for all data
    sq_diff   <- (X - mu[k])^2
    val       <- t(post.resp[, k]) %*% sq_diff
    sigma2[k] <- val / sum(post.resp[, k])
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
x <- seq(from = min(X)-1, to = max(X)+1, by = 0.1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0.0, .20), main=paste("Density plot of",K,"GMMs"), xlab = "x")
mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dnorm(x, mean=mu[k],sd=sqrt(sigma2[k]))
  mixture <- mixture + density * pi[k]
}
lines(x,mixture,col="red",lwd=2)