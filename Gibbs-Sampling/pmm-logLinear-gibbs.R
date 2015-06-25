#' Perform Gibbs sampling algorithm for Poisson Mixture Models (PMMs) used for 
#' modelling NGS RNA-Seq data. 
#' 
#' The parameter values are initialized using 'kmeans' algorithm. 
#' 
pmm.LL.gibbs <-  function(X, K=2, N.Sims=10000, burnin=5000, params){
  
  # Unwrap parameters from the 'params' list object
  conds     <- params$conds
  libSize   <- params$libSize
  libType   <- params$libType
  rm(params)
  
  N             <- NROW(X)                    # Length of the dataset
  q             <- NCOL(X)                    # Number of variables
  w             <- rowSums(X)                 # Overall expression for each object
  d             <- length(unique(conds))      # Total number of conditions
  r             <- as.vector(table(conds))    # Number of replicates in each condition
  
  post.resp     <- matrix(0, nrow=N, ncol=K)  # Posterior responsibilities
  pdf.w         <- matrix(0, nrow=N, ncol=K)  # PDF of each point on each cluster k
  lambdas       <- matrix(0, nrow=d, ncol=K)  # Matrix for holding estimated lambdas
  total.l       <- matrix(0, nrow=d, ncol=K)  # Store the sum of posterior means
  mean.mat      <- vector("list", K)          # List for holding the mean matrices l
  
  C.n           <- matrix(0, nrow=N, ncol=K)  # Mixture components
  C.matrix      <- matrix(0, nrow=N, ncol=K)  # Total Mixture components
  NLL           <- vector(mode="numeric")     # Hold NLL for all MCMC iterations
  x.k.bar       <- vector(length=K)           # Sample mean for each cluster k 
  lambda.draws  <- list()                     # Mean vector of each Poisson
  pi.draws      <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mixing Proportions
  #if (stephens) # Use Stephens algorithm for relabelling MCMC outputs
  #  postRespArr   <- array(0, dim=c(N.Sims-burnin, N, K)) # Post resp for each MCMC run
  
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
  
  ##===============================================
  # If 'params' not defined initialize parameters #
  ##===============================================
  # Initialize parameters using 'kmeans'
  initParams <- kmeansInit(X=X, 
                           K=K, 
                           w=w, 
                           s.dot=s.dot, 
                           conds=conds, 
                           lambdas=lambdas, 
                           eqProp=eqProp)
  Poisson       <- list()
  Poisson$l     <- initParams$lambdas           # Poisson mean vector for each cluster
  Poisson$Gamma <- list(shape.0=1, rate.0=1)    # Initialize Gamma hyperparameters
  pi.cur        <- initParams$pi.c              # Mixing proportions
  dir.a         <- rep(1/K, K)                  # Dirichlet concentration parameter

  for (t in 1:N.Sims){
    # Compute mean matrix using the estimated lambdas, normalization factors s and
    # the overall expression levels for each object w.
    for (k in 1:K){
      lambda.mat    <- matrix(rep(rep(Poisson$l[,k], times=r), each=N), nrow=N, ncol=q)
      mean.mat[[k]] <- w.mat * s.mat * lambda.mat
    }
    # Compute responsibilities
    res         <- compute.resp(X, pdf.w, K, Poisson, mean.mat, pi.cur)
    post.resp   <- res$post.resp
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, K, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k)
    # Update posterior mean
    Poisson$l   <- lambda.update(X, K, C.n, Poisson, w, s.dot)
    
    # Keep only the simulations after the burned in period has passed
    if (t > burnin){
      total.l                     <- total.l + Poisson$l
      NLL                         <- c(NLL, res$NLL)
      C.matrix                    <- C.matrix + C.n
      pi.draws[t - burnin,]       <- pi.cur
      lambda.draws[[t - burnin]]  <- Poisson$l
    }
  }
  
  # Object to keep input data
  dat         <- NULL
  dat$X       <- X
  dat$K       <- K
  dat$N       <- N
  dat$N.Sims  <- N.Sims
  dat$burnin  <- burnin
  
  # Object to hold all the MCMC draws
  draws       <- NULL
  draws$pi    <- pi.draws
  draws$l     <- lambda.draws
  
  # Object to hold the summaries for the parameters
  summary     <- NULL
  summary$pi  <- apply(pi.draws, 2, mean)     # Expected value of mix. prop.
  summary$l   <- total.l / ((N.Sims-burnin))
  # Add names to the estimated variables for clarity
  names(summary$pi)   <- paste("Clust", 1:K)
  colnames(summary$l) <- paste("Clust", 1:K)
  rownames(summary$l) <- conds.names
  
  #summary$l   <- apply(lambda.draws, 2, mean) # Expected value of mean
  summary$C   <- C.matrix / (N.Sims-burnin)   # Convert C.matrix to probs
  summary$NLL <- NLL
  
  # Object to hold the credible intervals for the parameters
  cred.interv     <- NULL
  cred.interv$pi  <- apply(pi.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  #cred.interv$l   <- apply(lambda.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  
  return(list(dat=dat, draws=draws, summary=summary, cred.interv=cred.interv))
}


# Compute the responsibilities
compute.resp <- function(X, pdf.w, K, Poisson, mean.mat, pi.cur){
  for (k in 1:K){
    pdf.w[,k] <- log(pi.cur[k]) + rowSums(dpois(X, lambda=mean.mat[[k]], log=TRUE))
  }
  # Calculate probabilities using the logSumExp trick for numerical stability
  Z           <- apply(pdf.w, 1, logSumExp)
  post.resp   <- pdf.w - Z
  post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
  NLL         <- -sum(Z)                    # Evaluate the NLL
  
  return(list(post.resp=post.resp, NLL=NLL))
}

# Update the mixture components 
c.n.update <- function(N, K, post.resp){
  c.i.draw <- matrix(0, nrow=N, ncol=K)
  for (i in 1:N){ # Sample one point from a multinomial i.e. ~ Discrete
    c.i.draw[i,] = rmultinom(1, 1, post.resp[i,])
    #c.i.draw[i] <- sample(1:clusters, size=1, prob=post.resp[i,],replace=TRUE)
  }
  return(c.i.draw)
}

# Update the mixing proportions
pi.update <- function(dir.a, N.k){
  a.n <- dir.a + N.k
  return(as.vector(rdirichlet(n=1, alpha=a.n))) # Sample from Dirichlet
}

# Update the posterior mean
lambda.update <- function(X, K, C.n, Poisson, w, s.dot){
  d           <- NROW(Poisson$l)            # Number of conditions
  lambda.post <- matrix(0, nrow=d, ncol=K)  # Matrix for holding estimated lambdas
  rate.up     <- colSums(C.n * w)           # Calculate RHS of denominator
  for (j in 1:d){
    beta.n    <- s.dot[j] * rate.up + Poisson$Gamma$rate.0
    X.j.      <- rowSums(as.matrix(X[,which(conds == (unique(conds))[j])]))
    alpha.n   <- colSums(C.n * matrix(rep(X.j., K), ncol=K)) + Poisson$Gamma$shape.0
    for (k in 1:K){
      lambda.post[j, k] <- rgamma(1, shape=alpha.n[k], rate=beta.n[k])
    }
  }
  return(lambda.post)
}
