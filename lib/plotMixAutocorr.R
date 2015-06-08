plotMixAutoc <- function(mcmc.obj, nams){
  ## --- Plot autocorrelation functions
  K     <- mcmc.obj$dat$K           # Number of clusters
  pars  <- length(mcmc.obj$draws)   # Number of parameters
  if (missing(nams)){
    nams  <- names(mcmc.obj$draws)  # Store names for showing in plots
  }
  par(mfrow=c(pars, K))             # Set plot frames
  
  for (i in 1:pars){            # Plot each parameter for ...
    for (j in 1:K) {            # each cluster k
      acf(mcmc.obj$draws[[i]][,j], main=paste(nams[i],j,sep="_"))
    }
  }
}