##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
source('../readData.R')
source('gmm1D-gibbs.R')

##===========================
# Initialize main variables #
##===========================
K           <- 3      # Number of clusters
N           <- 500    # Number of objects
N.Sims      <- 10000  # Set the number of simulations
burnin      <- 1000   # Set how many samples should be burned in
Normal      <- list() # Create a Normal object

##====================
# Generate the data  #
##====================
X <- gen.gaussian(N=N, K=K, pi=c(.4,.3,.3), mus=c(10,20,30), stds=c(2,3,1))

##=========================
# Initialize parameters   #
##=========================
cl              <- kmeans(X, K, nstart = 25) # Use Kmeans with random starts
C.n             <- cl$cluster   # get the mixture components
pi.cur          <- as.vector(table(C.n)/length(X)) # mixing proportions
dir.a           <- rep(1/K, K)  # Dirichlet concentration parameter vector
Normal$mu       <- as.vector(cl$centers) # Normal mean for each cluster
for (k in 1:K){
  Normal$Tau[k] <- 1/var(X[C.n==k]) # Normal precision for each cluster
}
Normal$Norm     <- list(mu.0 = 0, tau.0=1/100) # Normal hyperparameters
Normal$Gamma    <- list(a=1, b=1) # Gamma hyperparameters

##===================================
# Do inference using Gibbs sampling #
##===================================
gibbs <- gmm1D.gibbs(X, K, N.Sims, burnin, Normal, pi.cur, dir.a)

##=====================================
# Plot the data points and their pdfs #
##=====================================
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 0.1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0.0, .20), main=paste("Density plot of",K,"Gaussian Mixtures"),xlab="x")

mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dnorm(x, mean=mean(gibbs$mu.draws[,k]),sd=sqrt(1/mean(gibbs$tau.draws[,k])))
  mixture <- mixture + density * mean(gibbs$pi.draws[,k])
}
lines(x,mixture,col="red",lwd=2)