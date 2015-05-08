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
N       <- 500      # Total number of objects to create
K       <- 3        # Number of clusters
r1      <- rbinom(n=N, size=40, prob=.9)  # Total number of trials
r       <- matrix(r1, ncol=1)
X       <- gen.binomial(K=3, pi.c=c(.4,.3,.3), p=c(0.2,0.7, 0.5), r=r)

##=========================
# Initialize variables    #
##=========================
N       <- length(X)                  # Length of the dataset
cl      <- kmeans(X/r, K, nstart=25)  # Use Kmeans with random starts
C.n     <- cl$cluster                 # get the mixture components
p       <- as.vector(cl$centers)      # mean for each cluster
pi.c    <- as.vector(table(C.n)/length(X)) # mixing proportions

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
    pdf.w[,k] <- pi.c[k] * dbinom(X, size=r, prob=p[k])
  }
  post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalization
  
  ##========
  # M-Step #
  ##========
  for (k in 1:K){
    N.k       <- sum(post.resp[,k])     # Sum of responsibilities for cluster k
    pi.c[k]   <- N.k / N                # Update mixing proportions for cluster k
    p[k]      <- t(post.resp[,k]) %*% (X/r) / N.k # Update probabilities
  }
  
  # Evaluate the log likelihood
  # ln p(X|mu,S,p) = Sum_{n=1}^{N}(ln(Sum_{k=1}^{K}(p_k * N(x_n|mu_k, S_k))))
  logLik   <- sum(log(colSums(pdf.w)))
  if (abs(logLik-prevLogLik) < epsilon){ # Check for convergence.
    break
  }
  print(logLik)
} #End of Expectation Maximization loop.