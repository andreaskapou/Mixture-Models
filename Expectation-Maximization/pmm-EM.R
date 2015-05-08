##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory

##====================
# Generate the data  #
##====================
epsilon <- 1e-06    # Convergence paramater
K       <- 3        # Number of clusters
X       <- gen.poisson(K=K, pi.c=c(.3,.2,.5), lambdas=c(5,20,50))

##=========================
# Initialize variables    #
##=========================
N       <- length(X)                  # Length of the dataset
cl      <- kmeans(X, K, nstart = 25)  # Use Kmeans with random starts
C.n     <- cl$cluster                 # get the mixture components
pi.c    <- as.vector(table(C.n)/length(X)) # mixing proportions
lambdas <- as.vector(cl$centers)      # means for each Gaussian

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
    pdf.w[,k] <- pi.c[k] * dpois(X, lambda=lambdas[k])
  }
  post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalization
  
  ##========
  # M-Step #
  ##========
  for (k in 1:K){
    N.k         <- sum(post.resp[,k])     # Sum of responsibilities for cluster k
    pi.c[k]     <- N.k / N                # Update mixing proportions for cluster k
    lambdas[k]  <- (t(post.resp[, k]) %*% X) / N.k # Update mean for cluster k
  }
  # Check for convergence.
  logLik   <- sum(log(colSums(pdf.w)))
  if (abs(logLik-prevLogLik) < epsilon){
    break
  }
  print(logLik)
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