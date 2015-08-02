pmmMV.gibbs <-  function(X, K=2, N.Sims=10000, burnin=5000, params){
  ##================================================================
  # Function which uses Gibbs sampling so as to find the posterior #
  # of the Hierarchical Dirichlet Finite Mixture Model with        #
  # Gamma prior on lambda, using one dimensional dataset.          #
  ##================================================================
  
  N             <- NROW(X)                    # Length of the dataset
  D             <- NCOL(X)                    # Dimension of the dataset
  post.resp     <- matrix(0, nrow=N, ncol=K)  # Posterior responsibilities
  pdf.w         <- matrix(0, nrow=N, ncol=K)  # PDF of each point on each cluster k
  C.n           <- matrix(0, nrow=N, ncol=K)  # Mixture components
  C.matrix      <- matrix(0, nrow=N, ncol=K)  # Total Mixture components
  NLL           <- vector(mode="numeric")     # Hold NLL for all MCMC iterations
  x.k.bar       <- matrix(0, nrow=K, ncol=D)  # Sample mean vector for each cluster k
  total.l       <- matrix(0, nrow=K, ncol=D)  # Store the sum of posterior means
  lambda.draws  <- list()                     # Mean vector of each Poisson
  pi.draws      <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mixing Proportions
  
  ##===============================================
  # If 'params' not defined initialize parameters #
  ##===============================================
  if (missing(params)){
    Poisson       <- list()
    cl            <- kmeans(X, K, nstart = 25)    # Use Kmeans with random starts
    C             <- cl$cluster                   # Get the mixture components
    pi.cur        <- as.vector(table(C)/NROW(X))  # Mixing proportions
    dir.a         <- rep(1/K, K)                  # Dirichlet concentration parameter
    Poisson$l     <- cl$centers                   # Poisson mean for each cluster
    Poisson$Gamma <- list(shape.0=1, rate.0=1)    # Initialize Gamma hyperparameters
  }else{
    Poisson     <- params$Poisson
    pi.cur      <- params$pi.cur
    dir.a       <- params$dir.a
  }
  
  for (t in 1:N.Sims){
    # Compute responsibilities
    res         <- compute.resp(X, pdf.w, K, Poisson, pi.cur)
    post.resp   <- res$post.resp
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, K, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k) 
    for (k in 1:K){
      if (N.k[k] == 0){ # If there are no objects in the cluster
        x.k.bar[k,]  <- 0
      }else{ # Else, calculate sample mean for each cluster
        x.k.bar[k,] <- colMeans(X[C.n[,k] == 1,])
      }
    }
    # Update posterior mean
    Poisson$l   <- lambda.update(K, D, Poisson, N.k, x.k.bar)
    
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
  summary$l   <- total.l / ((N.Sims-burnin))  # Expected value of each mean vector
  # Add names to the estimated variables for clarity
  names(summary$pi)   <- paste("Clust", 1:K)
  rownames(summary$l) <- paste("Clust", 1:K)
  
  summary$C   <- C.matrix / (N.Sims-burnin)   # Convert C.matrix to probs
  summary$NLL <- NLL
  
  # Object to hold the credible intervals for the parameters
  cred.interv     <- NULL
  cred.interv$pi  <- apply(pi.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  
  return(list(dat=dat, draws=draws, summary=summary, cred.interv=cred.interv))
}


# Compute the responsibilities
compute.resp <- function(X, pdf.w, K, Poisson, pi.cur){
  N   <- NROW(X)  # Length of the dataset
  D   <- NCOL(X)  # Dimension of the dataset
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- log(pi.cur[k]) + rowSums(dpois(X, lambda=matrix(rep(Poisson$l[k,], each=N), ncol=D), log=TRUE))
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
lambda.update <- function(K, D, Poisson, N.k, x.k.bar){
  lambda.post <- matrix(0, nrow=K, ncol=D)  # Matrix for holding estimated lambdas
  for (k in 1:K){
    alpha.n   <- Poisson$Gamma$shape.0 + N.k[k] * x.k.bar[k,]
    beta.n    <- rep(Poisson$Gamma$rate.0 + N.k[k], D)
    lambda.post[k,] <- rgamma(D, shape=alpha.n, rate=beta.n)
  }
  return(lambda.post)
}
