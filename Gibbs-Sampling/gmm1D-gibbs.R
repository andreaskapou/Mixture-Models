gmm1D.gibbs <-  function(X, K=2, N.Sims, burnin, Normal, pi.cur, dir.a, logl=TRUE){
  ##===================================================================
  # Function which uses Gibbs sampling so as to find the posterior    #
  # of the Hierarchical Dirichlet Finite Mixture Model with           #
  # Normal-Gamma priors on mu and Tau, using one dimensional dataset  #
  ##===================================================================

  N             <- length(X)        # Length of the dataset
  post.resp     <- matrix(, N, K)   # Posterior responsibilities
  pdf.w         <- matrix(, N, K)   # PDF of each point on each cluster k
  C.n           <- matrix(, N, K)   # Mixture components
  x.k.bar       <- vector(length=K) # Sample mean for each cluster k 
  ssd.k         <- vector(length=K) # Sum of square differences on each cluster k
  mu.draws      <- matrix(, N.Sims-burnin, K) # Mean of each Gaussian
  tau.draws     <- matrix(, N.Sims-burnin, K) # Precision of each Gaussian
  pi.draws      <- matrix(, N.Sims-burnin, K) # Mixing Proportions

  for (t in 1:N.Sims){
    # Compute responsibilities
    post.resp   <- compute.resp(X, pdf.w, K, Normal, pi.cur, logl) 
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k) 
    for (k in 1:K){ # 
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
    
    # Keep only the simulations after the burn-in period has passed
    if (t > burnin){
      pi.draws[t - burnin,]     <- pi.cur
      mu.draws[t - burnin,]     <- Normal$mu
      tau.draws[t - burnin,]    <- Normal$Tau
    }
  }
  return(list(pi.draws=pi.draws, mu.draws=mu.draws, tau.draws=tau.draws))
}
  
# Compute the responsibilities
compute.resp <- function(X, pdf.w, K, Normal, pi.cur, logl){
  if (logl){
    for (k in 1:K) # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- log(pi.cur[k]) + dnorm(X, mean=Normal$mu[k], 
                                          sd=1/sqrt(Normal$Tau[k]), log=TRUE)
    post.resp   <-  pdf.w - apply(pdf.w,1,logSumExp) # Normalize the log probability
    post.resp   <- apply(post.resp, 2, exp) # Eponentiate to get actual probabilities
  }else{
    for (k in 1:K) # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- pi.cur[k] * dnorm(X, mean=Normal$mu[k], sd=1/sqrt(Normal$Tau[k]))
    post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  }
  return(post.resp)
}
# Update the mixture components 
c.n.update <- function(N, post.resp){
  c.i.draw <- matrix(, N, K)
  for (i in 1:N){ # Sample one point from a multinomial i.e. ~ Discrete
    c.i.draw[i,] = rmultinom(1, 1, post.resp[i,])
    #c.i.draw[i] <- which(as.vector(rmultinom(1,1,post.resp[i,])) == TRUE)
  }
  return(c.i.draw)
}
# Update the mixing proportions
pi.update <- function(dir.a, N.k){
  a.n <- dir.a + N.k
  return(as.vector(rdirichlet(n=1, alpha=a.n))) # Sample from Dirichlet
}
# Update the posterior precision
tau.update <- function(K, Normal, N.k, x.k.bar, ssd.k){
  tau.posterior <- vector(length=K)
  for (k in 1:K){
    alpha.n <- Normal$Gamma$a + N.k[k]/2
    beta.n  <- Normal$Gamma$b + ssd.k[k]/2 + ( 
          (N.k[k] * Normal$Norm$tau.0 * ((x.k.bar[k] - Normal$Norm$mu.0)^2)) / 
                                           (2 * (N.k[k] + Normal$Norm$tau.0)) )
    tau.posterior[k] <- rgamma(1, alpha.n, beta.n) # Sample from Gamma
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