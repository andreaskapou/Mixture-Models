##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
library(R.utils)
require(mvtnorm)
source('gmmMV-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)

##===========================
# Initialize main variables #
##===========================
K         <- 3      # Number of clusters
N         <- 200    # Number of objects
N.Sims    <- 10000  # Set the number of simulations
burnin    <- 1000   # Set how many samples should be burned in
Normal    <- list() # Create a Normal object

##====================
# Generate the data  #
##====================
#X        <- as.matrix(faithful)
mu        <- matrix(c(1.5,4,6,6,4,1.5),nrow=3,ncol=2) 
S         <- list()
S[[1]]    <- matrix(c(1,0,0,1),nrow=2,ncol=2)
S[[2]]    <- matrix(c(1,.7,.7,1),nrow=2,ncol=2)
S[[3]]    <- matrix(c(1,.5,.5,1),nrow=2,ncol=2)
X         <- gen.MV.gaussian(N=N, K=K, pi.c=c(.4,.3,.3), mus=mu, S=S)
D         <- NCOL(X) # Number of columns (variables - dimensions)
##=========================
# Initialize parameters   #
##=========================
cl        <- kmeans(X, K, nstart = 25)    # Use Kmeans with random starts
C.n       <- cl$cluster                   # get the mixture components
pi.cur    <- as.vector(table(C.n)/NROW(X))# Mixing proportions
dir.a     <- rep(1/K, K)                  # Dirichlet concentration parameter
logl      <- TRUE                         # If we want to compute log likel

invisible(readline(prompt="Run Gibbs sampling with diagonal covariance matrix..."))
for (k in 1:K){
  Normal$mu[[k]]  <- cl$centers[k,]       # Normal mean for each cluster
                                          # Normal precision matrix for each cluster
  Normal$Tau[[k]] <- 1/(apply(X[C.n==k,], MARGIN=2, FUN='sd')^2)
}
Normal$Norm   <- list(mu.0=colMeans(X), tau.0=rep(1/100,D)) # Normal hyperparameters 
Normal$Gamma  <- list(a=rep(1,D), b=apply(X, 2, sd)^2)      # Gamma hyperparameters

##===================================
# Do inference using Gibbs sampling #
##===================================
gibbs <- gmmMV.gibbs(X, K, N.Sims, burnin, Normal, pi.cur, dir.a, logl)

##=====================================
# Plot the data points and their pdfs #
##=====================================
invisible(readline(prompt="Press [enter] to show the plot"))
# Points to estimate the density and plot the contours
xpts <- seq(from=min(X[,1])-1,to=max(X[,1])+1,length.out=100)
ypts <- seq(from=min(X[,2])-1,to=max(X[,2])+1,length.out=100)
# Estimate the mixture contours
mixture.contour <- outer(xpts,ypts,function(x,y) {
  mixture <- matrix(data=0, length(xpts)*length(ypts), 1)
  for(k in 1:K){
    tau.k <- gibbs$total.Tau[[k]]/(N.Sims-burnin)
    mu.k  <- gibbs$total.mu[k,]/(N.Sims-burnin)
    mixture <- mixture + gibbs$pi.draws[k] * 
      dmvnorm(cbind(x,y),mean=mu.k,sigma=solve(tau.k))
  }
  return(mixture)
})
contour(xpts,ypts,mixture.contour,nlevels=10,drawlabel=FALSE,col="red",
        xlab="Eruption time",ylab="Waiting time",
        main=paste("Waiting vs Eruption time") )
points(X)


#####################################################################
#####################################################################


invisible(readline(prompt="Run Gibbs sampling with full covariance matrix..."))
##=========================
# Initialize parameters   #
##=========================
Normal    <- list() # Create a Normal object
for (k in 1:K){
  Normal$mu[[k]]  <- cl$centers[k,]       # Normal mean for each cluster
                                          # Normal precision matrix for each cluster
  Normal$Tau[[k]] <- solve(cov(X[C.n==k,]))
}
Normal$Norm    <- list(mu.0=colMeans(X), tau.0=1/100*diag(D)) 
Normal$Wishart <- list(v=D+2, k.0=0.01)

##===================================
# Do inference using Gibbs sampling #
##===================================
gibbs <- gmmMV.gibbs(X, K, N.Sims, burnin, Normal, pi.cur, dir.a, logl)

##=====================================
# Plot the data points and their pdfs #
##=====================================
invisible(readline(prompt="Press [enter] to show the plot"))
# Points to estimate the density and plot the contours
xpts <- seq(from=min(X[,1])-1,to=max(X[,1])+1,length.out=100)
ypts <- seq(from=min(X[,2])-1,to=max(X[,2])+1,length.out=100)

# Estimate the mixture contours
mixture.contour <- outer(xpts,ypts,function(x,y) {
  mixture <- matrix(data=0, length(xpts)*length(ypts), 1)
  for(k in 1:K){
    tau.k <- gibbs$total.Tau[[k]]/(N.Sims-burnin)
    mu.k  <- gibbs$total.mu[k,]/(N.Sims-burnin)
    mixture <- mixture + gibbs$pi.draws[k] * 
      dmvnorm(cbind(x,y),mean=mu.k,sigma=solve(tau.k))
  }
  return(mixture)
})
contour(xpts,ypts,mixture.contour,nlevels=10,drawlabel=FALSE,col="red",
        xlab="Eruption time",ylab="Waiting time",
        main=paste("Waiting vs Eruption time") )
points(X)