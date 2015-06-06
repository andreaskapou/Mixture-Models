#' Performs EM algorithm for univariate Gaussian Mixture Models (GMMs).
#' 
#' Notation and algorithm follows Bishop's book Ch.9, "Pattern Recognition 
#' and Machina Learning". BUT, it computes the Negative Log Likelihood (NLL):
#' - ln p(X|mu,S,p) = - Sum_{N}(ln(Sum_{K}(p_k * N(x_n|mu_k, S_k))))
#' 
#' This method works in different modes, depending on the parameters given. 
#' If no 'theta' parameter is given, then it initializes the parameters using
#' 'kmeans' algorithm. 

gmm.EM <- function(X, K=2, theta, epsilon=1e-10, maxIter=1000, isLog=TRUE, isDebug=FALSE){
  
  post.resp <- matrix(0, nrow=N, ncol=K)        # Hold responsibilities
  pdf.w     <- matrix(0, nrow=N, ncol=K)        # Hold weighted PDFs
  all.NLL   <- vector(mode="numeric")           # Hold NLL for all EM iterations
  NLL       <- 1e+40                            # Initialize Negative Log Likelihood
  
  # If 'theta' parameter is empty, we initialize parameters using 'kmeans'
  if (missing(theta)){
    N       <- length(X)                        # Length of the dataset
    cl      <- kmeans(X, K, nstart = 25)        # Use Kmeans with random starts
    C.n     <- cl$cluster                       # Get the mixture components
    mu      <- as.vector(cl$centers)            # Mean for each cluster
    pi.c    <- as.vector(table(C.n)/length(X))  # Mixing proportions
    s2      <- vector(length=K)                 # Variance for each cluster
    for (k in 1:K){
      s2[k] <- var(X[C.n==k])
    }
  }else{
    mu      <- theta$mu                         # Mean for each cluster
    s2      <- theta$s2                         # Variance for each cluster
    pi.c    <- theta$pi.c                       # Mixing proportions
  }
  
  if (isDebug){
    cat("Initial values:\n")
    cat("Mean:", mu, "\n")
    cat("Variance:", s2, "\n")
    cat("Mixing proportions:", pi.c, "\n")
    
    cat("Initial NLL:", NLL, "\n")
  }
  
  ##=========================================
  # Run Expectation Maximization  algorithm #
  ##=========================================
  for (i in 1:maxIter){                         # Loop until convergence
    prevNLL  <- NLL                             # Store NLL to check for convergence
    
    ##===================
    #       E-Step      #
    ##===================
    if (!isLog){
      # Calculate weighted PDF of each cluster for each data point
      for (k in 1:K){ 
        pdf.w[,k] <- pi.c[k] * dnorm(X, mean=mu[k], sd=sqrt(sigma2[k]))
      }
      Z           <- rowSums(pdf.w)             # Normalization constant
      post.resp   <- pdf.w / Z                  # Get responsibilites by normalization
      NLL         <- -sum(log(Z))               # Evaluate the NLL
    }else{
      for (k in 1:K){
        pdf.w[,k] <- log(pi.c[k]) + dnorm(X, mean=mu[k], sd=sqrt(sigma2[k]), log=TRUE)
      }
      # Calculate probabilities using the logSumExp trick for numerical stability
      Z           <- apply(pdf.w, 1, logSumExp)
      post.resp   <- pdf.w - Z
      post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
      NLL         <- -sum(Z)                    # Evaluate the NLL
    }
    
    ##===================
    #       M-Step      #
    ##===================
    N.k   <- colSums(post.resp)                 # Sum of responsibilities for each cluster
    pi.c  <- N.k/N                              # Update mixing proportions for each cluster
    mu    <- (t(X) %*% post.resp)/N.k           # Update mean for each cluster
    s2    <- (t(X^2) %*% post.resp)/N.k - mu^2  # Update variance for each cluster
    
    
    if (isDebug){
      cat("i:", i, "\n")
      cat("NLL:", NLL, "\n")
    }
    
    NLL.Diff  <- prevNLL - NLL                  # Compute NLL difference after ith iteration
    if (NLL.Diff < 0){
      stop("Negative log likelihood increases - Something is wrong!")
    }
    
    all.NLL   <- c(all.NLL, NLL)                # Keep all NLL in a vector  
    if (NLL.Diff < epsilon){                    # Check for convergence.
      break
    }
    
  } #End of Expectation Maximization loop.
  
  message("Total iterations: ", i, "\n")
  if (i == maxIter){
    message("Warning: EM did not converge with given maximum iterations!\n\n")
  }
  
  return(list(mu=mu, s2=s2, pi.c=pi.c, NLL=NLL, post.resp=post.resp, all.NLL=all.NLL))
}
