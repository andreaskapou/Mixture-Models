gmm1D.gibbs <-  function(X, K=2, N.Sims=10000, burnin=5000, params, logl=TRUE){
  ##===================================================================
  # Function which uses Gibbs sampling so as to find the posterior    #
  # of the Hierarchical Dirichlet Finite Mixture Model with           #
  # Normal-Gamma priors on mu and Tau, using one dimensional dataset  #
  ##===================================================================

  N             <- length(X)                  # Length of the dataset
  post.resp     <- matrix(0, nrow=N, ncol=K)  # Posterior responsibilities
  pdf.w         <- matrix(0, nrow=N, ncol=K)  # PDF of each point on each cluster k
  C.n           <- matrix(0, nrow=N, ncol=K)  # Mixture components
  C.matrix      <- matrix(0, nrow=N, ncol=K)  # Total Mixture components
  NLL           <- vector(mode="numeric")     # Hold NLL for all MCMC iterations
  x.k.bar       <- vector(length=K)           # Sample mean for each cluster k 
  ssd.k         <- vector(length=K)           # Sum of square differences on each cluster k
  mu.draws      <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mean of each Gaussian
  tau.draws     <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Precision of each Gaussian
  pi.draws      <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mixing Proportions
  
  ##===============================================
  # If 'Normal' not defined initialize parameters #
  ##===============================================
  if (missing(params)){
    Normal          <- list()                       # Create Normal object
    cl              <- kmeans(X, K, nstart = 25)    # Use Kmeans with random starts
    C               <- cl$cluster                   # Get the mixture components
    pi.cur          <- as.vector(table(C)/NROW(X))  # Mixing proportions
    dir.a           <- rep(1/K, K)                  # Dirichlet concentration parameter
    Normal$mu       <- as.vector(cl$centers)        # Normal mean for each cluster
    for (k in 1:K){
      Normal$Tau[k] <- 1/var(X[C==k])               # Normal precision for each cluster
    }
    Normal$Norm     <- list(mu.0=mean(X), tau.0=1/sd(X))  # Normal hyperparameters
    Normal$Gamma    <- list(shape.0=1, rate.0=1)          # Gamma hyperparameters
    
    ##=============================================================
    # Ordering constraint:                                        #
    #   Order all the values according to the mean parameter 'mu' #
    ##=============================================================
    Order       <- order(Normal$mu)
    pi.cur      <- pi.cur[Order]
    Normal$mu   <- Normal$mu[Order]
    Normal$Tau  <- Normal$Tau[Order]
  } else{
    Normal          <- params$Normal
    pi.cur          <- params$pi.cur
    dir.a           <- params$dir.a
  }

  for (t in 1:N.Sims){  # Start Gibbs sampling
    # Compute responsibilities
    res         <- compute.resp(X, pdf.w, K, Normal, pi.cur, logl)
    post.resp   <- res$post.resp
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, K, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k)
    for (k in 1:K){
      if (N.k[k] == 0){ # If there are no objects in the cluster
        x.k.bar[k]  <- 0
        ssd.k[k]    <- 0
      }else{ # Else, calculate sample mean and variance for each cluster
        x.k.bar[k]  <- mean(X[C.n[,k] == 1])
        ssd.k[k]    <- sum((X[C.n[,k] == 1] - x.k.bar[k])^2)
      }
    }
    # Update posterior precision
    Normal$Tau <- tau.update(K, Normal, N.k, x.k.bar, ssd.k)
    # Update posterior mean
    Normal$mu  <- mu.update(K, Normal, N.k, x.k.bar)
    
    ##=============================================================
    # Ordering constraint:                                        #
    #   Order all the values according to the mean parameter 'mu' #
    ##=============================================================
    Order       <- order(Normal$mu)
    C.n         <- C.n[, Order]
    pi.cur      <- pi.cur[Order]
    Normal$mu   <- Normal$mu[Order]
    Normal$Tau  <- Normal$Tau[Order]
    
    # Keep only the simulations after the burn-in period has passed
    if (t > burnin){
      NLL                       <- c(NLL, res$NLL)
      C.matrix                  <- C.matrix + C.n
      pi.draws[t - burnin,]     <- pi.cur
      mu.draws[t - burnin,]     <- Normal$mu
      tau.draws[t - burnin,]    <- Normal$Tau
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
  draws$mu    <- mu.draws
  draws$tau   <- tau.draws
  
  # Object to hold the summaries for the parameters
  summary     <- NULL
  summary$pi  <- apply(pi.draws, 2, mean)   # Expected value of mix. prop.
  summary$mu  <- apply(mu.draws, 2, mean)   # Expected value of mean
  summary$tau <- apply(tau.draws, 2, mean)  # Expected value of variance
  summary$C   <- C.matrix / (N.Sims-burnin) # Convert C.matrix to probs
  summary$NLL <- NLL
  
  # Object to hold the credible intervals for the parameters
  cred.interv     <- NULL
  cred.interv$pi  <- apply(pi.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  cred.interv$mu  <- apply(mu.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  cred.interv$tau <- apply(tau.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  
  return(list(dat=dat, draws=draws, summary=summary, cred.interv=cred.interv))
}


# Compute the responsibilities
compute.resp <- function(X, pdf.w, K, Normal, pi.cur, logl){
  if (logl){
    for (k in 1:K){ # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- log(pi.cur[k]) + dnorm(X, mean=Normal$mu[k], 
                                          sd=1/sqrt(Normal$Tau[k]), log=TRUE)
    }
    # Calculate probabilities using the logSumExp trick for numerical stability
    Z           <- apply(pdf.w, 1, logSumExp)
    post.resp   <- pdf.w - Z
    post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
    NLL         <- -sum(Z)                    # Evaluate the NLL
  }else{
    for (k in 1:K){ # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- pi.cur[k] * dnorm(X, mean=Normal$mu[k], sd=1/sqrt(Normal$Tau[k]))
    }
    Z           <- rowSums(pdf.w)             # Normalization constant
    post.resp   <- pdf.w / Z                  # Get responsibilites by normalization
    NLL         <- -sum(log(Z))               # Evaluate the NLL
  }
  return(list(post.resp=post.resp, NLL=NLL))
}

# Update the mixture components 
c.n.update <- function(N, K, post.resp){
  c.i.draw <- matrix(0, nrow=N, ncol=K)
  for (i in 1:N){ # Sample one point from a multinomial i.e. ~ Discrete
    c.i.draw[i,] = rmultinom(n=1, size=1, post.resp[i,])
    #c.i.draw[i] <- which(as.vector(rmultinom(1,1,post.resp[i,])) == TRUE)
  }
  return(c.i.draw)
}

# Update the mixing proportions
pi.update <- function(dir.a, N.k){
  a.n <- dir.a + N.k
  return(as.vector(MCMCpack::rdirichlet(n=1, alpha=a.n))) # Sample from Dirichlet
}

# Update the posterior precision
tau.update <- function(K, Normal, N.k, x.k.bar, ssd.k){
  tau.posterior <- vector(length=K)
  for (k in 1:K){
    alpha.n <- Normal$Gamma$shape.0 + N.k[k]/2
    beta.n  <- Normal$Gamma$rate.0 + ssd.k[k]/2 + ( 
          (N.k[k] * Normal$Norm$tau.0 * ((x.k.bar[k] - Normal$Norm$mu.0)^2)) / 
                                           (2 * (N.k[k] + Normal$Norm$tau.0)) )
    tau.posterior[k] <- rgamma(1, shape=alpha.n, rate=beta.n) # Sample from Gamma
  }
  return(tau.posterior)
}

# Update the posterior mean
mu.update <- function(K, Normal, N.k, x.k.bar){
  mu.posterior <- vector(length=K)
  for (k in 1:K){
    mu.n    <- (Normal$Norm$mu.0 * Normal$Norm$tau.0 + N.k[k] * x.k.bar[k]) / 
                  (Normal$Norm$tau.0 + N.k[k])
    tau.n <- Normal$Tau[k] * (Normal$Norm$tau.0 + N.k[k])
    mu.posterior[k] <- rnorm(1, mean=mu.n, sd=1/sqrt(tau.n)) # Sample from Normal
  }
  return(mu.posterior)
}
