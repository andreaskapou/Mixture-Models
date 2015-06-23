#' Performs EM algorithm for multivariate Gaussian Mixture Models (GMMs).
#' 
#' Notation and algorithm follows Bishop's book Ch.9, "Pattern Recognition 
#' and Machina Learning". BUT, it computes the Negative Log Likelihood (NLL):
#' - ln p(X|mu,S,p) = - Sum_{N}(ln(Sum_{K}(p_k * N(x_n|mu_k, S_k))))
#' 
#' This method works in different modes, depending on the parameters given. 
#' If no 'theta' parameter is given, then it initializes the parameters using
#' 'kmeans' algorithm. 
##==================================================================================
gmmMV.EM <- function(X, K=2, theta, epsilon=1e-10, maxIter=1000, isLog=TRUE, isDebug=FALSE){
  
  N         <- NROW(X)                          # Length of the dataset
  D         <- NCOL(X)                          # Dimension of the dataset
  post.resp <- matrix(0, nrow=N, ncol=K)        # Hold responsibilities
  pdf.w     <- matrix(0, nrow=N, ncol=K)        # Hold weighted PDFs
  all.NLL   <- vector(mode="numeric")           # Hold NLL for all EM iterations
  NLL       <- 1e+40                            # Initialize Negative Log Likelihood
  
  # If 'theta' parameter is empty, we initialize parameters using 'kmeans'
  if (missing(theta)){
    cl      <- kmeans(X, K, nstart = 25)        # Use Kmeans with random starts
    C.n     <- cl$cluster                       # Get the mixture components
    mu      <- cl$centers                       # Mean for each cluster
    pi.c    <- as.vector(table(C.n)/NROW(X))    # Mixing proportions
    Sigma   <- list(length=K)                   # Variance for each cluster
    for (k in 1:K){
      Sigma[[k]]  <- cov(X[C.n==k,])
    }
  }else{
    mu      <- theta$mu                         # Mean for each cluster
    Sigma   <- theta$Sigma                      # Variance for each cluster
    pi.c    <- theta$pi.c                       # Mixing proportions
  }
  
  if (isDebug){
    cat("Initial values:\n")
    cat("Mean:", mu, "\n")
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
        pdf.w[,k] <- pi.c[k] * mvtnorm::dmvnorm(X, mean=mu[k,], sigma=Sigma[[k]], log=F)
      }
      Z           <- rowSums(pdf.w)             # Normalization constant
      post.resp   <- pdf.w / Z                  # Get responsibilites by normalization
      NLL         <- -sum(log(Z))               # Evaluate the NLL
    }else{
      for (k in 1:K){
        pdf.w[,k] <- log(pi.c[k]) + mvtnorm::dmvnorm(X, mean=mu[k,], sigma=Sigma[[k]], log=T)
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
    N.k           <- colSums(post.resp)         # Sum of responsibilities for each cluster
    pi.c          <- N.k / N                    # Update mixing proportions for each cluster
    for (k in 1:K){
      mu[k,]      <- (t(post.resp[,k]) %*% X) / N.k[k]            # Update mu
      X.cen       <- sweep(X, MARGIN=2, STATS=mu[k,], FUN="-")
      Sigma[[k]]  <- t(X.cen) %*% (X.cen*post.resp[,k]) / N.k[k]  # Update Sigma
    }
    
    
    NLL.Diff  <- prevNLL - NLL                  # Compute NLL difference after ith iteration
    if (NLL.Diff < 0){
      message("Negative log likelihood increases - Something is wrong!\n")
      message("Finishing EM...!")
      break
    }
    if (isDebug){
      cat("i:", i, "\t")
      cat("NLL:", NLL, "\t\t")
      cat("NLL-diff:", NLL.Diff, "\n")
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
  
  # Add names to the estimated variables for clarity
  names(pi.c)   <- paste("Clust", 1:K)
  rownames(mu)  <- paste("Clust", 1:K)
  names(Sigma)  <- paste("Clust", 1:K)
  
  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility.
  labels <- unlist(apply(post.resp, 1, function(x) which(x == max(x, na.rm = TRUE))[1]))
  
  ##===========================
  # Perform model selection   #
  ##===========================
  numParams <- (K-1) + K*D + K*(D*(D+1)/2)  # Total number of parameters i.e. pi.c + mu + S
  
  BIC <- 2*NLL + numParams*log(N) # BIC = -2*ln(L) + params*ln(N)
  AIC <- 2*NLL + 2*numParams      # AIC = -2*ln(L) + 2*params
  
  entropy <- -sum(post.resp * log(post.resp), na.rm=TRUE)
  ICL <- BIC + entropy            # Integrated Complete Likelihood criterion
  
  ##=========================
  # Create a GMM-MV object  #
  ##=========================
  results           <- list()
  results$X         <- X
  results$K         <- K
  results$N         <- N
  results$D         <- D
  results$postResp  <- post.resp
  results$labels    <- labels
  results$pi.c      <- pi.c
  results$mu        <- mu
  results$Sigma     <- Sigma
  results$NLL       <- NLL
  results$all.NLL   <- all.NLL
  results$BIC       <- BIC
  results$AIC       <- AIC
  results$ICL       <- ICL
  
  class(results) <- "GMM-MV"
  return(results)
}
