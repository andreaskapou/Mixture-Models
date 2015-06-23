#' Performs EM algorithm for Binomial distributed Probit Regression mixture model.
#' 
#' We compute the Negative Log Likelihood (NLL):
#' - ln p(X|theta) = - Sum_{N}(ln(Sum_{K}(p_k * BinProb(x_n|theta))))
##=================================================================================

binProbReg.EM <- function(X, K=2, params, epsilon=1e-4, maxIter=1000, isDebug=FALSE){
  
  N         <- length(X)                        # Length of the dataset
  post.resp <- matrix(0, nrow=N, ncol=K)        # Hold responsibilities
  pdf.w     <- matrix(0, nrow=N, ncol=K)        # Hold weighted PDFs
  all.NLL   <- vector(mode="numeric")           # Hold NLL for all EM iterations
  NLL       <- 1e+40                            # Initialize Negative Log Likelihood
  
  pi.c      <- params$pi.c                      # Mixing proportions
  theta     <- params$theta                     # Polynomial coefficients
  
  if (isDebug){
    cat("Initial values:\n")
    cat("Theta:", theta, "\n")
    cat("Mixing proportions:", pi.c, "\n")
    
    cat("Initial NLL:", NLL, "\n")
  }
  
  ##=========================================
  # Run Expectation Maximization  algorithm #
  ##=========================================
  for (t in 1:maxIter){                         # Loop until convergence
    prevNLL  <- NLL                             # Store NLL to check for convergence
    
    ##===================
    #       E-Step      #
    ##===================
    # Calculate weighted PDF of each cluster for each data point
    for (k in 1:K){
      for (i in 1:N){
        pdf.w[i,k] <- log(pi.c[k]) + binProbRegLik(theta=theta[,k], D=X[[i]], mode=1)
      }
    }
    # Calculate probabilities using the logSumExp trick for numerical stability
    Z           <- apply(pdf.w, 1, logSumExp)
    post.resp   <- pdf.w - Z
    post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
    NLL         <- -sum(Z)                    # Evaluate the NLL
    
    ##===================
    #       M-Step      #
    ##===================
    N.k   <- colSums(post.resp)                 # Sum of responsibilities for each cluster
    pi.c  <- N.k/N                              # Update mixing proportions for each cluster
    for (k in 1:K){
      # Update the parameters a, b, c of the polynomial using Conjugate-Gradient method
      theta[,k] <- optim(par=theta[,k],
                         fn=sumBinPRL,
                         gr=sumDerBinPRL, X, post.resp[,k],
                         method="CG",
                         control = list(fnscale=-1, maxit = 20) )$par
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
  
  message("Total iterations: ", t, "\n")
  if (t == maxIter){
    message("Warning: EM did not converge with given maximum iterations!\n\n")
  }
  
  # Add names to the estimated variables for clarity
  names(pi.c)     <- paste("Clust", 1:K)
  colnames(theta) <- paste("Clust", 1:K)
  
  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility.
  labels <- unlist(apply(post.resp, 1, function(x) which(x == max(x, na.rm = TRUE))[1]))
  
  ##===========================
  # Perform model selection   #
  ##===========================
  numParams <- (K-1) + K*NROW(theta)  # Total number of parameters i.e. pi.c + theta
  
  BIC <- 2*NLL + numParams*log(N) # BIC = -2*ln(L) + params*ln(N)
  AIC <- 2*NLL + 2*numParams      # AIC = -2*ln(L) + 2*params
  
  entropy <- -sum(post.resp * log(post.resp), na.rm=TRUE)
  ICL <- BIC + entropy            # Integrated Complete Likelihood criterion
  
  ##======================
  # Create a BPR object  #
  ##======================
  results           <- list()
  results$X         <- X
  results$K         <- K
  results$N         <- N
  results$postResp  <- post.resp
  results$labels    <- labels
  results$pi.c      <- pi.c
  results$theta     <- theta
  results$NLL       <- NLL
  results$all.NLL   <- all.NLL
  results$BIC       <- BIC
  results$AIC       <- AIC
  results$ICL       <- ICL
  
  class(results) <- "BPR"
  return(results)
}
