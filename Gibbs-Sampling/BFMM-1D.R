##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
source('../readData.R')
source('bmm1D-gibbs.R')
source('logSumExp.R')

##===========================
# Initialize main variables #
##===========================
K           <- 3      # Number of clusters
N           <- 500    # Number of objects
N.Sims      <- 10000  # Set the number of simulations
burnin      <- 1000   # Set how many samples should be burned in
Binom       <- list() # Create a Binomial object

##====================
# Generate the data  #
##====================
r <- rbinom(n=N, size=50, prob=.8)
X <- gen.binomial(N=N, K=K, pi=c(.4,.3,.3), p=c(0.2,0.9, 0.4), r=r)

##=========================
# Initialize variables    #
##=========================
cl          <- kmeans(X/r, K, nstart = 25)  # Use Kmeans with random starts
C.n         <- cl$cluster                   # Get the mixture components
pi.cur      <- as.vector(table(C.n)/length(X)) # Mixing proportions
dir.a       <- rep(1/K, K)                  # Dirichlet concentration parameter
Binom$p     <- as.vector(cl$centers)        # Binomial probability for each cluster
Binom$Beta  <- list(a=1, b=1)               # Initialize Beta hyperparameters
Binom$r     <- matrix(r, ncol=1)            # Total number of trials for each object

##===================================
# Do inference using Gibbs sampling #
##===================================
gibbs <- bmm1D.gibbs(X, K, N.Sims, burnin, Binom, pi.cur, dir.a)

##=====================================
# Plot the data points and their pdfs #
##=====================================
# Create x points from min(X) to max(X)
x <- seq(from = min(X)-1, to = max(X)+1, by = 1)
hist(X, breaks = 22, freq=FALSE, col="lightblue", xlim=c(min(X)-1,max(X)+1),
     ylim = c(0.0, .10), main=paste("Density plot of",K,"Binomial Mixtures"),xlab="x")

mixture = 0.0
for (k in 1:K){
  # Calculate the estimated density
  density <- dbinom(x, size=40, prob=mean(gibbs$p.draws[,k]))
  mixture <- mixture + density * mean(gibbs$pi.draws[,k])
}
lines(x,mixture,col="red",lwd=2)