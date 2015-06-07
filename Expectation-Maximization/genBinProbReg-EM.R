#' Performs EM algorithm for Binomial distributed Probit Regression mixture model.
#' 
#' We compute the Negative Log Likelihood (NLL):
#' - ln p(X|theta) = - Sum_{N}(ln(Sum_{K}(p_k * BinProb(x_n|theta))))

genBinProbReg.EM <- function(X, K=2, params, epsilon=1e-4, maxIter=1000, isDebug=FALSE){
  
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
        pdf.w[i,k] <- log(pi.c[k]) + genBinomProbRegrLik(theta=theta[,k], D=X[[i]], mode=1)
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
                         fn=sumGenBinomProbRegrLik,
                         gr=sumGenDerBinomProbRegrLik, X, post.resp[,k],
                         method="CG",
                         control = list(fnscale=-1, maxit = 5) )$par
    }
    
    
    if (isDebug){
      cat("i:", t, "\n")
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
  
  message("Total iterations: ", t, "\n")
  if (t == maxIter){
    message("Warning: EM did not converge with given maximum iterations!\n\n")
  }
  
  return(list(theta=theta, pi.c=pi.c, NLL=NLL, post.resp=post.resp, all.NLL=all.NLL))
}
