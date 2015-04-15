##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
library(R.utils)
source('pmm1D-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)

##===========================
# Initialize main variables #
##===========================
K           <- 3      # Number of clusters
N           <- 500    # Number of objects
N.Sims      <- 10000  # Set the number of simulations
burnin      <- 1000   # Set how many samples should be burned in
Poisson     <- list() # Create a Poisson object

##====================
# Generate the data  #
##====================
X <- gen.poisson(N=N, K=K, pi.c=c(.3,.2,.5), lambdas=c(7,26,60))

##=========================
# Initialize variables    #
##=========================
cl            <- kmeans(X, K, nstart = 25)    # Use Kmeans with random starts
C.n           <- cl$cluster                   # Get the mixture components
pi.cur        <- as.vector(table(C.n)/NROW(X))# Mixing proportions
dir.a         <- rep(1/K, K)                  # Dirichlet concentration parameter
Poisson$l     <- as.vector(cl$centers)        # Poisson mean for each cluster
Poisson$Gamma <- list(a=1, b=1)               # Initialize Gamma hyperparameters
logl          <- TRUE                         # If we want to compute log likel

##===================================
# Do inference using Gibbs sampling #
##===================================
gibbs <- pmm1D.gibbs(X, K, N.Sims, burnin, Poisson, pi.cur, dir.a, logl)

##=====================================
# Plot the data points and their pdfs #
##=====================================
invisible(readline(prompt="Press [enter] to show the plot"))
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0.0, .10), main=paste("Density plot of",K,"Poisson Mixtures"),xlab="x")

mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dpois(x, lambda=mean(gibbs$lambda.draws[,k]))
  mixture <- mixture + density * mean(gibbs$pi.draws[,k])
}
lines(x,mixture,col="red",lwd=2)