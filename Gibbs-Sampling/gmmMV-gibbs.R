gmmMV.gibbs <-  function(X, K=2, N.Sims, burnin, Normal, pi.cur, dir.a, logl){
  ##=====================================================================
  # Function which uses Gibbs sampling so as to find the posterior      #
  # of the Hierarchical Dirichlet Finite Mixture Model with             #
  # Normal-Gamma priors on mu and Tau, using multi-dimensional dataset  #
  ##=====================================================================
  
  N         <- NROW(X)  # Number of rows (observations)
  D         <- NCOL(X)  # Number of columns (variables - dimensions)  
  
  pdf.w     <- matrix(0, nrow=N, ncol=K)  # PDF of each point on each cluster k
  post.resp <- matrix(0, nrow=N, ncol=K)  # Matrix for the posterior responsibilites
  x.k.bar   <- matrix(0, nrow=K, ncol=D)  # Sample mean vector for each cluster k
  if (is.vector(Normal$Tau[[1]])){        # If diagonal covariance matrix
    ssd.k   <- matrix(0, nrow=K, ncol=D)  # Sum of square differences on cluster k
  }else{
    ssd.k   <- list()     # Sum of square differences on each cluster k
  }
  total.mu  <- matrix(0, nrow=K, ncol=D)  # Store the sum of posterior means
  total.Tau <- list()                     # Store the sum of posterior precisions
  for (k in 1:K){
    total.Tau[[k]] <- matrix(0, nrow=D, ncol=D)
  }
  pi.draws  <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mixing Proportions
  mu.draws  <- list()     # Mean vector of each Gaussian
  Tau.draws <- list()     # Precision matrix of each Gaussian
  
  
  for (t in 1:N.Sims){    # Start Gibbs sampling
    # Compute responsibilities
    post.resp   <- compute.resp(X, pdf.w, K, Normal, pi.cur, logl)
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, K, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k) 
    
    for (k in 1:K){ # Calculate sample mean and SSD for each cluster
      x.k.bar[k,] <- colMeans(X[C.n[,k] == 1,])
                    # apply(X[C.n[,k]==1,], MARGIN=2, FUN='sd')^2 * N.k[k]
      if (is.vector(Normal$Tau[[1]])){  # If diagonal covariance matrix
        ssd.k[k,]   <- rowSums(apply(X[C.n[,k] == 1,],1,'-',x.k.bar[k,])^2)
      }else{
        ssd.k[[k]] <- apply(X[C.n[,k] == 1,],1,'-',x.k.bar[k,]) %*% 
                              t(apply(X[C.n[,k] == 1,],1,'-',x.k.bar[k,]))
      }
    }
    if (is.vector(Normal$Tau[[1]])){  # If diagonal covariance matrix
      # Update posterior precision 
      Tau.cur <- tau.update.diag(K, D, Normal, N.k, x.k.bar, ssd.k)
      # Update posterior mean
      mu.cur  <- mu.update.diag(K, D, Normal, N.k, x.k.bar)
    }else{
      # Update posterior precision
      Tau.cur <- tau.update.full(K, D, Normal, N.k, x.k.bar, ssd.k)
      # Update posterior mean
      mu.cur  <- mu.update.full(K, D, Normal, N.k, x.k.bar)
    }
    
    # Keep only the simulations after the burn-in period has passed
    if (t > burnin){
      total.mu  <- total.mu + mu.cur
      if (is.vector(Normal$Tau[[1]])){
        for (k in 1:K){
          total.Tau[[k]]  <- total.Tau[[k]] + diag(Tau.cur[k,])
        }
      }else{
        for (k in 1:K){
          total.Tau[[k]]  <- total.Tau[[k]] + Tau.cur[[k]]
        }
      }
      mu.draws[[t - burnin]]  <- mu.cur
      Tau.draws[[t - burnin]] <- Tau.cur
      pi.draws[t - burnin,]   <- pi.cur
    }
  }
  return(list(pi.draws=pi.draws, mu.draws=mu.draws, Tau.draws=Tau.draws,
                                  total.mu=total.mu, total.Tau=total.Tau))
}

# Compute the responsibilities
compute.resp <- function(X, pdf.w, K, Normal, pi.cur, logl){
  if (logl){
    if (is.vector(Normal$Tau[[1]])){ # If diagonal covariance matrix
      for (k in 1:K)  # Calculate the PDF of each cluster for each data point
        pdf.w[,k] <- log(pi.cur[k]) + dmvnorm(X, mean=Normal$mu[[k]],
                                        sigma=diag(1/Normal$Tau[[k]]), log=TRUE)
    }else{
      for (k in 1:K)  # Calculate the PDF of each cluster for each data point
        pdf.w[,k] <- log(pi.cur[k]) + dmvnorm(X, mean=Normal$mu[[k]], 
                                        sigma=solve(Normal$Tau[[k]]), log=TRUE)
    }
    post.resp   <-  pdf.w - apply(pdf.w,1,logSumExp) # Normalize the log probability
    post.resp   <- apply(post.resp, 2, exp) # Eponentiate to get actual probabilities
  }else{
    if (is.vector(Normal$Tau[[1]])){ # If diagonal covariance matrix
      for (k in 1:K)  # Calculate the PDF of each cluster for each data point
        pdf.w[,k] <- pi.cur[k] * dmvnorm(X, mean=Normal$mu[[k]],
                                         sigma=diag(1/Normal$Tau[[k]]))
    }else{
      for (k in 1:K)  # Calculate the PDF of each cluster for each data point
        pdf.w[,k] <- pi.cur[k] * dmvnorm(X, mean=Normal$mu[[k]], 
                                         sigma=solve(Normal$Tau[[k]]))
    }
    post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion
  }
  return(post.resp)
}
# Update the mixture components
c.n.update <- function(N, K, post.resp){
  c.i.draw <- matrix(, nrow=N, ncol=K)
  for (i in 1:N)  # Sample one point from a multinomial i.e. ~ Discrete
    c.i.draw[i,] = rmultinom(1, 1, post.resp[i,])
  return(c.i.draw)
}
# Update the mixing proportions
pi.update <- function(dir.a, N.k){
  a.n <- dir.a + N.k
  return(as.vector(rdirichlet(n=1, alpha=a.n))) # Sample from Dirichlet
}
# Update the posterior precision diagonal covariance matrix
tau.update.diag <- function(K, D, Normal, N.k, x.k.bar, ssd.k){
  tau.posterior <- matrix(0, nrow=K, ncol=D)
  for (k in 1:K){
    alpha.n <- Normal$Gamma$a + N.k[k]/2
    beta.n  <- Normal$Gamma$b + ssd.k[k,]/2 + ( 
      (N.k[k] * Normal$Norm$tau.0 * ((x.k.bar[k,] - Normal$Norm$mu.0)^2)) / 
        (2 * (N.k[k] + Normal$Norm$tau.0)) )
    tau.posterior[k,] <- rgamma(D, alpha.n, beta.n) # Sample from Gamma
  }
  return(tau.posterior)
}
# Update the posterior precision full covariance matrix
tau.update.full <- function(K, D, Normal, N.k, x.k.bar, ssd.k){
  tau.posterior <- list()
  for (k in 1:K){
    v.n   <- Normal$Wishart$v + N.k[k]
    f     <- (x.k.bar[k,]-Normal$Norm$mu.0)%*%(t(x.k.bar[k,]-Normal$Norm$mu.0))
    Tau.n <- Normal$Norm$tau.0 + ssd.k[[k]] + ( (N.k[k]*Normal$Wishart$k.0*f) / 
                                                  (Normal$Wishart$k.0 + N.k[k]) )
    tau.posterior[[k]] <- rwish(v.n, solve(Tau.n)) # Sample from Wishart
  }
  return(tau.posterior)
}
# Update the posterior mean diagonal covariance matrix
mu.update.diag <- function(K, D, Normal, N.k, x.k.bar){
  mu.posterior <- matrix(0, nrow=K, ncol=D)
  for (k in 1:K){
    t.n.k   <- Normal$Norm$tau.0 + N.k[k]
    mu.n    <- (Normal$Norm$mu.0*Normal$Norm$tau.0 + N.k[k]*x.k.bar[k,]) / t.n.k
    Tau.n   <- Normal$Tau[[k]] * t.n.k
    mu.posterior[k,] <- rnorm(D, mean=mu.n, sd=1/sqrt(Tau.n)) # Sample from Normal
  }
  return(mu.posterior)
}
# Update the posterior mean full covariance matrix
mu.update.full <- function(K, D, Normal, N.k, x.k.bar){
  mu.posterior <- matrix(0, nrow=K, ncol=D)
  for (k in 1:K){
    t.n.k   <- Normal$Wishart$k.0 + N.k[k]
    mu.n    <- (Normal$Norm$mu.0*Normal$Wishart$k.0+N.k[k]*x.k.bar[k,]) / t.n.k
    Tau.n   <- Normal$Tau[[k]] * t.n.k
    mu.posterior[k,] <- rmvnorm(1,mean=mu.n,sigma=solve(Tau.n)) # Sample MV-Normal
  }
  return(mu.posterior)
}