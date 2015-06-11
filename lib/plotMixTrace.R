plotMixTrace <- function(mcmc.obj, nams){
  # Plot traces of sampled parameters during MCMC
  K     <- mcmc.obj$dat$K           # Number of clusters
  pars  <- length(mcmc.obj$draws)   # Number of parameters
  if (missing(nams)){
    nams  <- names(mcmc.obj$draws)  # Store names for showing in plots
  }
  par(mfrow=c(pars, K))             # Set plot frames
  
  for (i in 1:pars){            # Plot each parameter for ...
    for (j in 1:K) {            # each cluster k
      plot(mcmc.obj$draws[[i]][,j], main=paste(nams[i],j,sep="_"), type="l")
    }
  }
}
