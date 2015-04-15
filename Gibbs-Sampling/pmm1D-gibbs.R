pmm1D.gibbs <-  function(X, K=2, N.Sims, burnin, Poisson, pi.cur, dir.a, logl=TRUE){
  ##================================================================
  # Function which uses Gibbs sampling so as to find the posterior #
  # of the Hierarchical Dirichlet Finite Mixture Model with        #
  # Gamma prior on lambda, using one dimensional dataset.          #
  ##================================================================
  
  N             <- length(X)                  # Length of the dataset
  post.resp     <- matrix(0, nrow=N, ncol=K)  # Posterior responsibilities
  pdf.w         <- matrix(0, nrow=N, ncol=K)  # PDF of each point on each cluster k
  C.n           <- matrix(0, nrow=N, ncol=K)  # Mixture components
  x.k.bar       <- vector(length=K)           # Sample mean for each cluster k 
  lambda.draws  <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mean of each Poisson
  pi.draws      <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mixing Proportions
  
  for (t in 1:N.Sims){
    # Compute responsibilities
    post.resp   <- compute.resp(X, pdf.w, K, Poisson, pi.cur, logl)
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, K, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k) 
    for (k in 1:K){
      if (N.k[k] == 0){ # If there are no objects in the cluster
        x.k.bar[k]  <- 0
      }else{ # Else, calculate sample mean for each cluster
        x.k.bar[k]  <- mean(X[C.n[,k] == 1])
      }
    }
    # Update posterior mean
    Poisson$l   <- lambda.update(K, Poisson, N.k, x.k.bar)

    # Keep only the simulations after the burned in period has passed
    if (t > burnin){
      pi.draws[t - burnin,]     <- pi.cur
      lambda.draws[t - burnin,] <- Poisson$l
    }
  }
  return(list(pi.draws=pi.draws, lambda.draws=lambda.draws))
}

# Compute the responsibilities
compute.resp <- function(X, pdf.w, K, Poisson, pi.cur, logl){
  if (logl){
    for (k in 1:K) # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- log(pi.cur[k]) + dpois(X, lambda=Poisson$l[k], log=TRUE)
    post.resp   <-  pdf.w - apply(pdf.w,1,logSumExp) # Normalize the log probability
    post.resp   <- apply(post.resp, 2, exp) # Eponentiate to get actual probabilities
  }else{
    for (k in 1:K) # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- pi.cur[k] * dpois(X, lambda=Poisson$l[k])
    post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  }
  return(post.resp)
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
lambda.update <- function(K, Poisson, N.k, x.k.bar){
  lambda.posterior <- vector(length=K)
  for (k in 1:K){
    alpha.n   <- Poisson$Gamma$a + N.k[k] * x.k.bar[k]
    beta.n    <- Poisson$Gamma$b + N.k[k]
    lambda.posterior[k] <- rgamma(1, shape=alpha.n, rate=beta.n)
  }
  return(lambda.posterior)
}