bmmMV.gibbs <-  function(X, r, K=2, N.Sims=10000, burnin=5000, params){
  ##================================================================
  # Function which uses Gibbs sampling so as to find the posterior #
  # of the Hierarchical Dirichlet Finite Mixture Model with        #
  # Beta prior on probability , using one dimensional dataset.     #
  ##================================================================
  
  N             <- NROW(X)                    # Length of the dataset
  D             <- NCOL(X)                    # Dimension of the dataset
  post.resp     <- matrix(0, nrow=N, ncol=K)  # Posterior responsibilities
  pdf.w         <- matrix(0, nrow=N, ncol=K)  # PDF of each point on each cluster k
  C.n           <- matrix(0, nrow=N, ncol=K)  # Mixture components
  C.matrix      <- matrix(0, nrow=N, ncol=K)  # Total Mixture components
  NLL           <- vector(mode="numeric")     # Hold NLL for all MCMC iterations
  sum.x.k       <- matrix(0, nrow=K, ncol=D)  # Sample sum for each cluster k 
  diff.k        <- matrix(0, nrow=K, ncol=D)  # Sum of difference for each cluster k 
  total.p       <- matrix(0, nrow=K, ncol=D)  # Store the sum of posterior successes
  p.draws       <- list()                     # Success prob of each trial
  pi.draws      <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mixing Proportions
  
  ##===============================================
  # If 'params' not defined initialize parameters #
  ##===============================================
  if (missing(params)){
    Binom       <- list()
    cl          <- kmeans(X/r, K, nstart = 25)  # Use Kmeans with random starts
    C           <- cl$cluster                   # Get the mixture components
    pi.cur      <- as.vector(table(C)/NROW(X))  # Mixing proportions
    dir.a       <- rep(1/K, K)                  # Dirichlet concentration parameter
    Binom$p     <- cl$centers                   # Binomial probability for each cluster
    Binom$Beta  <- list(a=1, b=1)               # Initialize Beta hyperparameters
  }else{
    Binom       <- params$Binom
    pi.cur      <- params$pi.cur
    dir.a       <- params$dir.a
  }
  
  for (t in 1:N.Sims){  # Start Gibbs sampling
    # Compute responsibilities
    res         <- compute.resp(X, r, pdf.w, K, Binom, pi.cur, logl)
    post.resp   <- res$post.resp
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, K, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k)
    for (k in 1:K){ 
      if (N.k[k] == 0){ # If there are no objects in the cluster
        sum.x.k[k,] <- 0
        diff.k[k,]  <- 0
      }else{ # Else, calculate sum and sum of difference for each cluster
        sum.x.k[k,] <- colSums(X[C.n[,k] == 1,])
        diff.k[k,]  <- colSums(r[C.n[,k] == 1,] - X[C.n[,k] == 1,])
      }
    }
    # Update posterior mean
    Binom$p <- p.update(K, D, Binom, sum.x.k, diff.k)
    
    # Keep only the simulations after the burned in period has passed
    if (t > burnin){
      total.p                 <- total.p + Binom$p
      NLL                     <- c(NLL, res$NLL)
      C.matrix                <- C.matrix + C.n
      pi.draws[t - burnin,]   <- pi.cur
      p.draws[[t - burnin]]   <- Binom$p
    }
  }
  
  # Object to keep input data
  dat         <- NULL
  dat$X       <- X
  dat$r       <- r
  dat$K       <- K
  dat$N       <- N
  dat$N.Sims  <- N.Sims
  dat$burnin  <- burnin
  
  # Object to hold all the MCMC draws
  draws       <- NULL
  draws$pi    <- pi.draws
  draws$p     <- p.draws
  
  # Object to hold the summaries for the parameters
  summary     <- NULL
  summary$pi  <- apply(pi.draws, 2, mean)   # Expected value of mix. prop.
  summary$p   <- total.p / ((N.Sims-burnin))  # Expected value of each mean vector
  # Add names to the estimated variables for clarity
  names(summary$pi)   <- paste("Clust", 1:K)
  rownames(summary$p) <- paste("Clust", 1:K)
  
  summary$C   <- C.matrix / (N.Sims-burnin) # Convert C.matrix to probs
  summary$NLL <- NLL
  
  # Object to hold the credible intervals for the parameters
  cred.interv     <- NULL
  cred.interv$pi  <- apply(pi.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  
  return(list(dat=dat, draws=draws, summary=summary, cred.interv=cred.interv))
}

# Compute the responsibilities
compute.resp <- function(X, r, pdf.w, K, Binom, pi.cur, logl){
  N   <- NROW(X)  # Length of the dataset
  D   <- NCOL(X)  # Dimension of the dataset
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- log(pi.cur[k]) + rowSums(dbinom(X, size=r, prob=matrix(rep(Binom$p[k,],each=N),ncol=D), log=TRUE))
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
  }
  return(c.i.draw)
}

# Update the mixing proportions
pi.update <- function(dir.a, N.k){
  a.n <- dir.a + N.k
  return(as.vector(rdirichlet(n=1, alpha=a.n))) # Sample from Dirichlet
}

# Update the posterior mean
p.update <- function(K, D, Binom, sum.x.k, diff.k){
  p.post <- matrix(0, nrow=K, ncol=D)  # Matrix for holding estimated p
  for (k in 1:K){
    alpha.n   <- Binom$Beta$a + sum.x.k[k,]
    beta.n    <- Binom$Beta$b + diff.k[k,]
    p.post[k,] <- rbeta(D, shape1=alpha.n, shape2=beta.n)
  }
  return(p.post)
}