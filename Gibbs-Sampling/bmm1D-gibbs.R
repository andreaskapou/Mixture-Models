bmm1D.gibbs <-  function(X, K=2, N.Sims, burnin, Binom, pi.cur, dir.a){
  ##================================================================
  # Function which uses Gibbs sampling so as to find the posterior #
  # of the Hierarchical Dirichlet Finite Mixture Model with        #
  # Beta prior on probability , using one dimensional dataset.     #
  ##================================================================
  
  N             <- length(X)        # Length of the dataset
  post.resp     <- matrix(, N, K)   # Posterior responsibilities
  pdf.w         <- matrix(, N, K)   # PDF of each point on each cluster k
  C.n           <- matrix(, N, K)   # Mixture components
  sum.x.k       <- vector(length=K) # Sample sum for each cluster k 
  diff.k        <- vector(length=K) # Sum of difference for each cluster k 
  p.draws       <- matrix(, N.Sims-burnin, K) # Success probability of each Binomial
  pi.draws      <- matrix(, N.Sims-burnin, K) # Mixing Proportions
  
  for (t in 1:N.Sims){
    # Compute responsibilities
    post.resp   <- compute.resp(X, pdf.w, K, Binom, pi.cur)
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k)
    for (k in 1:K){ 
      if (N.k[k] == 0){ # If there are no objects in the cluster
        sum.x.k[k]  <- 0
        diff.k[k]   <- 0
      }else{ # Else, calculate sum and sum of difference for each cluster
        sum.x.k[k]  <- sum(X[C.n[,k] == 1])
        diff.k[k]   <- sum((Binom$r[C.n[,k] == 1] - X[C.n[,k] == 1]))
      }
    }
    # Update posterior mean
    Binom$p <- p.update(K, Binom, sum.x.k, diff.k)
    
    # Keep only the simulations after the burned in period has passed
    if (t > burnin){
      pi.draws[t - burnin,]   <- pi.cur
      p.draws[t - burnin,]    <- Binom$p
    }
  }
  return(list(pi.draws=pi.draws, p.draws=p.draws))
}

# Compute the responsibilities
compute.resp <- function(X, pdf.w, K, Binom, pi.cur){
  for (k in 1:K) # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- pi.cur[k] * dbinom(X, size=Binom$r, prob=Binom$p[k])
  post.resp <- (pdf.w / rowSums(pdf.w)) # Get responsibilites by normalizarion
  return(post.resp)
}
# Update the mixture components 
c.n.update <- function(N, post.resp){
  c.i.draw <- matrix(, N, K)
  for (i in 1:N){ # Sample one point from a multinomial i.e. ~ Discrete
    c.i.draw[i,] = rmultinom(1, 1, post.resp[i,])
  }
  return(c.i.draw)
}
# Update the mixing proportions
pi.update <- function(dir.a, N.k){
  a.n <- dir.a + N.k
  return(as.vector(rdirichlet(n=1, alpha=a.n))) # Sample from Dirichlet
}
# Update the posterior mean
p.update <- function(K, Binom, sum.x.k, diff.k){
  p.posterior <- vector(length=K)
  for (k in 1:K){
    alpha.n   <- Binom$Beta$a + sum.x.k[k]
    beta.n    <- Binom$Beta$b + diff.k[k]
    p.posterior[k] <- rbeta(1, shape1=alpha.n, shape2=beta.n)
  }
  return(p.posterior)
}