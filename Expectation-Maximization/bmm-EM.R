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
N       <- 500        # Total number of objects to create
K       <- 3          # Number of clusters
r1      <- rbinom(n=N, size=40, prob=.8)  # Total number of trials
r       <- matrix(r1, ncol=1)
X       <- gen.binomial(K=3, pi=c(.4,.3,.3), p=c(0.2,0.7, 0.5), r=r)

##=========================
# Initialize variables    #
##=========================
N       <- length(X)              # Length of the dataset
cl      <- kmeans(X/r, K, nstart = 25) # Use Kmeans with random starts
C.n     <- cl$cluster             # get the mixture components
p       <- as.vector(cl$centers)  # mean for each cluster
pi      <- as.vector(table(C.n)/length(X)) # mixing proportions

##===============================
# Run Expectation Maximization  #
##===============================
# Hold responsibilities that component k takes on explaining the observation x.n
post.resp   <- matrix(, N, K)
log.likel   <- 0

for (i in 1:1000){  # Loop until convergence
  ##===============================================
  ## Expectation Step
  pdf.w       <- matrix(, N, K) # Hold the PDF of each point on each cluster k
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- pi[k] * dbinom(X, size=r, prob=p[k])
  }
  post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  
  ##===============================================
  ## Maximization Step
  prev.log.likel <- log.likel # Store to check for convergence
  for (k in 1:K){
    # Update mixing proportions
    pi[k]     <- sum(post.resp[, k] / sum(post.resp))
    # Update probabilities by taking the weighted average of *all* data points.
    val       <- t(post.resp[, k]) %*% (X/r)
    p[k]      <- val / sum(post.resp[,k])
  }
  # Check for convergence.
  log.likel   <- sum(log(sum(pdf.w)))
  if (abs(log.likel-prev.log.likel) < epsilon){
    break
  }
  print(log.likel)
} #End of Expectation Maximization loop.