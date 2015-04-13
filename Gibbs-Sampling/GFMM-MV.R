##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
require(mvtnorm)
source('../readData.R')
source('gmmMV-gibbs.R')
source('logSumExp.R')

##===========================
# Initialize main variables #
##===========================
K         <- 2      # Number of clusters
N.Sims    <- 10000  # Set the number of simulations
burnin    <- 1000   # Set how many samples should be burned in
Normal    <- list() # Create a Normal object

##====================
# Generate the data  #
##====================
X   <- as.matrix(faithful)
D   <- NCOL(X)      # Number of columns (variables - dimensions) 
##=========================
# Initialize parameters   #
##=========================
cl        <- kmeans(X, K, nstart = 25)    # Use Kmeans with random starts
C.n       <- cl$cluster                   # get the mixture components
pi.cur    <- as.vector(table(C.n)/NROW(X))# Mixing proportions
dir.a     <- rep(1/K, K)                  # Dirichlet concentration parameter
logl      <- TRUE                         # If we want to compute log likel
for (k in 1:K){
  Normal$mu[[k]]  <- cl$centers[k,]       # Normal mean for each cluster
                                          # Normal precision matrix for each cluster
  Normal$Tau[[k]] <- 1/(apply(X[C.n==k,], MARGIN=2, FUN='sd')^2)
}

#Normal$Norm    <- list(mu.0=colMeans(X), tau.0=1/100*diag(D)) 
Normal$Gamma   <- list(a=rep(1,D), b=apply(X, 2, sd)^2) 
Normal$Norm     <- list(mu.0=rep(0,D), tau.0=rep(1/100,D))  # Normal hyperparameters
#Normal$Gamma    <- list(a=rep(1,D), b=rep(1,D))             # Gamma hyperparameters

##===================================
# Do inference using Gibbs sampling #
##===================================
invisible(readline(prompt="Run Gibbs sampling with diagonal covariance matrix..."))
gibbs <- gmmMV.gibbs(X, K, N.Sims, burnin, Normal, pi.cur, dir.a, logl)

invisible(readline(prompt="Press [enter] to show the plot"))
##=====================================
# Plot the data points and their pdfs #
##=====================================
# Points to estimate the density and plot the contours
xpts <- seq(from=min(X[,1])-1,to=max(X[,1])+1,length.out=100)
ypts <- seq(from=min(X[,2])-10,to=max(X[,2])+10,length.out=100)
# Estimate the mixture contours
mixture.contour <- outer(xpts,ypts,function(x,y) {
  mixture <- matrix(data=0, length(xpts)*length(ypts), 1)
  for(k in 1:K){
    tau.k <- gibbs$total.Tau[[k]]/(N.Sims-burnin)
    mu.k  <- gibbs$total.mu[k,]/(N.Sims-burnin)
    #print(mu.k)
    mixture <- mixture + gibbs$pi.draws[k] * 
      dmvnorm(cbind(x,y),mean=mu.k,sigma=solve(tau.k))
  }
  return(mixture)
})

contour(xpts,ypts,mixture.contour,nlevels=10,drawlabel=FALSE,col="red",
        xlab="Eruption time",ylab="Waiting time",
        main=paste("Waiting vs Eruption time - Old Faithful") )
points(faithful)

#####################################################################
#####################################################################


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
invisible(readline(prompt="Run Gibbs sampling with full covariance matrix..."))
gibbs <- gmmMV.gibbs(X, K, N.Sims, burnin, Normal, pi.cur, dir.a, logl)

invisible(readline(prompt="Press [enter] to show the plot"))
##=====================================
# Plot the data points and their pdfs #
##=====================================
# Points to estimate the density and plot the contours
xpts <- seq(from=min(X[,1])-1,to=max(X[,1])+1,length.out=100)
ypts <- seq(from=min(X[,2])-10,to=max(X[,2])+10,length.out=100)
# Estimate the mixture contours
mixture.contour <- outer(xpts,ypts,function(x,y) {
  mixture <- matrix(data=0, length(xpts)*length(ypts), 1)
  for(k in 1:K){
    tau.k <- gibbs$total.Tau[[k]]/(N.Sims-burnin)
    mu.k  <- gibbs$total.mu[k,]/(N.Sims-burnin)
    #print(mu.k)
    mixture <- mixture + gibbs$pi.draws[k] * 
      dmvnorm(cbind(x,y),mean=mu.k,sigma=solve(tau.k))
  }
  return(mixture)
})

contour(xpts,ypts,mixture.contour,nlevels=10,drawlabel=FALSE,col="red",
        xlab="Eruption time",ylab="Waiting time",
        main=paste("Waiting vs Eruption time - Old Faithful") )
points(faithful)