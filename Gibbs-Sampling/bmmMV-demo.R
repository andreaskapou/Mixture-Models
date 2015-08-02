##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
library(coda)
library(label.switching)
library(R.utils)
source('bmmMV-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)
set.seed(12345)

##=============================================
# Initialize main variables and generate data #
##=============================================
N       <- 1000         # Number of data points
K       <- 2            # Number of clusters
r1      <- rbinom(n=N, size=40, prob=.9)  # Total number of trials 1st dim
r2      <- rbinom(n=N, size=60, prob=.9)  # Total number of trials 2nd dim
r       <- as.matrix(cbind(r1, r2))
X       <- gen.MV.binomial(N=N, r=r)

##=========================
# Initialize parameters   #
##=========================
N.Sims      <- 10000        # Set the number of simulations
burnin      <- 5000         # Set how many samples should be discarded
dir.a       <- rep(1, K)    # Dirichlet concentration parameter
pi.cur      <- rep(1/K, K)  # Initialize mixing proportions for each cluster

Binom       <- list()       # Create object of type Binomial
Binom$p     <- matrix(c(.1,.2,.3,.4), ncol=2) # Initialize probs of success for each cluster
Binom$Beta  <- list(a=1, b=1)     # Initialize Beta hyperparameters

params      <- list(Binom=Binom, pi.cur=pi.cur, dir.a=dir.a)

##=======================================================
# Do inference using Gibbs sampling, explicitly giving  #
# initial parameters.                                   #
##=======================================================
gibbs <- bmmMV.gibbs(X=X, 
                     r=r, 
                     K=K, 
                     N.Sims=N.Sims, 
                     burnin=burnin, 
                     params=params)

##=================================================================
# Do inference using Gibbs sampling, without initial parameters,  #
# the method will initilize parameters using k-means              #
##=================================================================
gibbs.kmeans <- bmmMV.gibbs(X=X,
                            r=r,
                            K=K, 
                            N.Sims=N.Sims, 
                            burnin=burnin)

