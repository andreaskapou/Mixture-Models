#' Performs EM algorithm for Poisson Mixture Models (PMMs) used for modelling 
#' NGS RNA-Seq data. 
#' 
#' The parameter values are initialized using 'kmeans' algorithm. 
#' 
pmm.LL.EM <- function(X, K=2, theta, epsilon=1e-06, maxIter=1000, isDebug=FALSE){
  
  # Unwrap parameters from the 'theta' list object
  conds     <- theta$conds
  eqProp    <- theta$eqProp
  libSize   <- theta$libSize
  libType   <- theta$libType
  rm(theta)
  
  N         <- NROW(X)                      # Length of the dataset
  q         <- NCOL(X)                      # Number of variables
  w         <- rowSums(X)                   # Overall expression for each object
  d         <- length(unique(conds))        # Total number of conditions
  r         <- as.vector(table(conds))      # Number of replicates in each condition
  
  post.resp <- matrix(0, nrow=N, ncol=K)    # Hold responsibilities
  pdf.w     <- matrix(0, nrow=N, ncol=K)    # Hold weighted PDFs
  lambdas   <- matrix(0, nrow=d, ncol=K)    # Matrix for holding estimated lambdas
  all.NLL   <- vector(mode="numeric")       # Hold NLL for all EM iterations
  mean.mat  <- vector("list", K)            # List for holding the mean matrices l
  NLL       <- 1e+40                        # Initialize Negative Log Likelihood
  
  # Grouping columns of X in order of condition (all replicates put together)
  o.ycols   <- order(conds)       # Order of conditions
  X         <- X[,o.ycols]        # Order the observations X accordingly
  conds     <- conds[o.ycols]     # Order the observation vector accordingly
  rm(o.ycols)
  
  # Make sure X is an N x q matrix and assign unique names to X and conds
  X <- as.matrix(X, nrow=N, ncol=q)
  
  if(length(rownames(X)) == 0){   # If matrix X has no row names
    rn <- 1:nrow(X)
  }else if(length(rownames(X)) > 0){ 
    rn <- rownames(X)
  }
  rownames(X) <- rn               # Assign names to each row of X
  conds.names <- unique(conds)    # Get unique condition names
  
  # Compute the library size normalization factors for each variable
  s <- normFactors(X, libSize=libSize, libType=libType)
  
  # Sum of s for all replicates l on each condition d
  s.dot <- rep(NA, d) 
  for (j in 1:d){
    s.dot[j] <- sum( s[which(conds == unique(conds)[j])] )
  }
  
  # Create matrices of dimension N x q, for faster computations
  w.mat     <- matrix(rep(w, times=q), nrow=N, ncol=q)
  s.mat     <- matrix(rep(s, each=N) , nrow=N, ncol=q)
  
  # Initialize parameters using 'kmeans'
  initParams <- kmeansInit(X=X, 
                           K=K, 
                           w=w, 
                           s.dot=s.dot, 
                           conds=conds, 
                           lambdas=lambdas, 
                           eqProp=eqProp)
  pi.c    <- initParams$pi.c
  lambdas <- initParams$lambdas
  
  if (isDebug){
    cat("Initial values:\n")
    cat("Mixing proportions:", pi.c, "\n")
    cat("Initial NLL:", NLL, "\n")
  }
  
  
  ##=========================================
  # Run Expectation Maximization  algorithm #
  ##=========================================
  for (t in 1:maxIter){                         # Loop until convergence
    prevNLL  <- NLL                             # Store NLL to check for convergence
    
    # Compute mean matrix using the estimated lambdas, normalization factors s and
    # the overall expression levels for each object w.
    for (k in 1:K){
      lambda.mat    <- matrix(rep(rep(lambdas[,k], times=r), each=N), nrow=N, ncol=q)
      mean.mat[[k]] <- w.mat * s.mat * lambda.mat
    }
    
    ##===================
    #       E-Step      #
    ##===================
    for (k in 1:K){
      pdf.w[,k] <- log(pi.c[k]) + rowSums(dpois(X, lambda=mean.mat[[k]], log=TRUE))
    }
    # Calculate probabilities using the logSumExp trick for numerical stability
    Z           <- apply(pdf.w, 1, logSumExp)
    post.resp   <- pdf.w - Z
    post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
    NLL         <- -sum(Z)                    # Evaluate the NLL
    
    
    ##===================
    #       M-Step      #
    ##===================
    if(eqProp){
      pi.c  <- rep(1/K, K)                    # Equal mixing proportions 
    }else{
      N.k     <- colSums(post.resp)           # Sum of responsibilities for each cluster
      pi.c    <- N.k/N                        # Update mixing proportions for each cluster
    }                                        
                                              # Update unknown parameters lambda
    lambdas <- calcLambdas(X=X,         
                           w=w, 
                           s.dot=s.dot, 
                           conds=conds, 
                           postResp=post.resp, 
                           lambdas=lambdas)
    
    
    NLL.Diff  <- prevNLL - NLL                # Compute NLL difference after ith iteration
    if (NLL.Diff < 0){
      message("Negative log likelihood increases - Something is wrong!\n")
      message("Finishing EM...!")
      break
    }
    if (isDebug){
      cat("t:", t, "\t")
      cat("NLL:", NLL, "\t\t")
      cat("NLL-diff:", NLL.Diff, "\n")
    }
    all.NLL   <- c(all.NLL, NLL)              # Keep all NLL in a vector  
    if (NLL.Diff < epsilon){                  # Check for convergence.
      break
    }
  } #End of Expectation Maximization loop.
  
  message("Total iterations: ", t, "\n")
  if (t == maxIter){
    message("Warning: EM did not converge with given maximum iterations!\n\n")
  }
  
  # Add names to the estimated variables for clarity
  names(pi.c)       <- paste("Clust", 1:K)
  colnames(lambdas) <- paste("Clust", 1:K)
  rownames(lambdas) <- conds.names
  
  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility.
  labels <- unlist(apply(post.resp, 1, function(x) which(x == max(x, na.rm = TRUE))[1]))
  
  ##===========================
  # Perform model selection   #
  ##===========================
  if (eqProp){
    numParams <- (K-1) + K*(d-1)  # Total number of parameters i.e. pi.c + lambda
  }else{
    numParams <- K*(d-1)          # Total number of parameters i.e. lambda
  }
  
  BIC <- 2*NLL + numParams*log(N) # BIC = -2*ln(L) + params*ln(N)
  AIC <- 2*NLL + 2*numParams      # AIC = -2*ln(L) + 2*params
  
  entropy <- -sum(post.resp * log(post.resp), na.rm=TRUE)
  ICL <- BIC + entropy            # Integrated Complete Likelihood criterion
  
  ##===============================
  # Create a PMMLogLinear object  #
  ##===============================
  results           <- list()
  results$X         <- X
  results$K         <- K
  results$w         <- w
  results$s         <- s
  results$N         <- N
  results$q         <- q
  results$d         <- d
  results$r         <- r
  results$postResp  <- post.resp
  results$labels    <- labels
  results$conds     <- conds
  results$libSize   <- libSize
  results$libType   <- libType
  results$pi.c      <- pi.c
  results$lambdas   <- lambdas
  results$NLL       <- NLL
  results$all.NLL   <- all.NLL
  results$BIC       <- BIC
  results$AIC       <- AIC
  results$ICL       <- ICL
  
  class(results) <- "PMMLogLinear"
  return(results)
}
