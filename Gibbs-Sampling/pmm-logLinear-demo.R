##=========================================================
# Script for running Poisson Log Linear Mixture Models    #
# using EM to fit the model to the data.                  #
##=========================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(MCMCpack)
library(coda)
library(label.switching)
library(R.utils)
library(edgeR)
source("pmm-logLinear-gibbs.R")
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)           # Set seed for reproducible results

##=============================================
# Generate the data and initialize variables  #
##=============================================
N         <- 1000               # Number of objects
sim       <- poisMixSim(N=N, libsize="A", separation="low")
X         <- sim$X              # Data points, Nxq matrix
conds     <- sim$conds          # Different biological conditions

K         <- 4                  # Number of clusters
N.Sims    <- 10000              # Set the number of simulations
burnin    <- 5000               # Set how many samples should be discarded
libSize   <- TRUE               # Take library size into consideration
libType   <- "TMM"              # Method to compute normalization factors
# Available: "TC", "MED", "DESEQ", "TMM"

# Wrap all the parameters in a list
params     <- list(conds=conds,
                   libSize=libSize,
                   libType=libType)

##=======================================================
# Run PMM-LL-Gibbs, parameters initialized by 'kmeans'  #
##=======================================================
fit.pmmLL.gibbs <- pmm.LL.gibbs(X=X,
                                K=K,
                                N.Sims=N.Sims,
                                burnin=burnin,
                                params=params)
