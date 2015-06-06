##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory

set.seed(1)
##====================
# Generate the data  #
##====================
epsilon <- 1e-06    # Convergence paramater
K       <- 3        # Number of clusters
X       <- gen.gaussian(K=K, pi.c=c(.4,.3,.3), mus=c(0,2,5), stds=c(1,1,1))

#gm<-normalmixEM(X,k=K,lambda=pi.c,mu=mu,sigma=sigma2)

##=========================
# Initialize variables    #
##=========================
N       <- length(X)                  # Length of the dataset
cl      <- kmeans(X, K, nstart = 25)  # Use Kmeans with random starts
C.n     <- cl$cluster                 # get the mixture components
mu      <- as.vector(cl$centers)      # means for each Gaussian
pi.c    <- as.vector(table(C.n)/length(X)) # mixing proportions
sigma2  <- vector(length=K)           # Variance for each Gaussian
for (k in 1:K){
  sigma2[k] <- var(X[C.n==k])
}
isLog       <- FALSE
post.resp   <- matrix(, N, K)         # Hold responsibilities
pdf.w       <- matrix(, N, K)         # Hold PDF of each point on each cluster k
logLik      <- 0                      # Initialize log likelihood

##===============================
# Run Expectation Maximization  #
##===============================
for (i in 1:100){                    # Loop until convergence
  prevLogLik  <- logLik               # Store to check for convergence
  
  ##========
  # E-Step #
  ##========
  if (!isLog){
    for (k in 1:K){ # Calculate the weighted PDF of each cluster for each data point
      pdf.w[,k] <- pi.c[k] * dnorm(X, mean=mu[k], sd=sqrt(sigma2[k]))
    }
    post.resp   <- pdf.w / rowSums(pdf.w)             # Get responsibilites by normalization
    
    # ln p(X|mu,S,p) = Sum_{n=1}^{N}(ln(Sum_{k=1}^{K}(p_k * N(x_n|mu_k, S_k))))
    logLik      <- sum(log(colSums(pdf.w)))           # Evaluate the log likelihood
  }else{
    for (k in 1:K){
      pdf.w[,k] <- log(pi.c[k]) + dnorm(X, mean=mu[k], sd=sqrt(sigma2[k]), log=TRUE)
    }
    post.resp   <- pdf.w - apply(pdf.w, 1, logSumExp) # Normalize the log probability
    post.resp   <- apply(post.resp, 2, exp)           # Exponentiate to get actual probabilities
    
    logLik      <- sum(log(colSums(exp(pdf.w))))      # Evaluate the log likelihood
  }
  
  ##========
  # M-Step #
  ##========
  N.k <- colSums(post.resp)                           # Sum of responsibilities for each cluster
  pi.c <- N.k / N                                     # Update mixing proportions for each cluster
  mu <- (t(X) %*% post.resp) / N.k                    # Update mean for each cluster
  sigma2 <- (t(X^2) %*% post.resp) / N.k - mu^2       # Update variance
  #for (k in 1:K){
  #  N.k       <- sum(post.resp[,k])                   
  #  pi.c[k]   <- N.k / N                              
  #  mu[k]     <- (post.resp[,k] %*% X) / N.k          
  #  sigma2[k] <- (post.resp[,k] %*% (X-mu[k])^2)/N.k
  #}
  
  if (abs(logLik - prevLogLik) < epsilon){            # Check for convergence.
    break
  }
  print(logLik)
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
  mixture <- mixture + density * pi.c[k]
}
lines(x,mixture,col="red",lwd=2)