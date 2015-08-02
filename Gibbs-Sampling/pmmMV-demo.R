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
source('pmmMV-gibbs.R')
sourceDirectory("../lib", modifiedOnly=FALSE)
set.seed(12345)

##=============================================
# Initialize main variables and generate data #
##=============================================
N       <- 1000         # Number of data points
K       <- 2            # Number of clusters
X       <- gen.MV.poisson(N=N)

##=========================
# Initialize variables    #
##=========================
N.Sims        <- 10000        # Set the number of simulations
burnin        <- 5000         # Set how many samples should be discarded
dir.a         <- rep(1, K)    # Dirichlet concentration parameter
pi.cur        <- rep(1/K, K)  # Initialize mixing proportions for each cluster

Poisson       <- list()       # Create an object of type Poisson
Poisson$l     <- matrix(c(1,2,3,4),ncol=2)     # Initialize means for each cluster
Poisson$Gamma <- list(shape.0=1, rate.0=0.01)  # Gamma hyperparameters

params        <- list(Poisson=Poisson, pi.cur=pi.cur, dir.a=dir.a)

##=======================================================
# Do inference using Gibbs sampling, explicitly giving  #
# initial parameters.                                   #
##=======================================================
gibbs <- pmmMV.gibbs(X=X, 
                     K=K, 
                     N.Sims=N.Sims, 
                     burnin=burnin, 
                     params=params)

##=================================================================
# Do inference using Gibbs sampling, without initial parameters,  #
# the method will initilize parameters using k-means              #
##=================================================================
gibbs.kmeans <- pmmMV.gibbs(X=X, 
                            K=K, 
                            N.Sims=N.Sims, 
                            burnin=burnin)

