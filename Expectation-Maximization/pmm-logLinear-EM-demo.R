##=========================================================
# Script for running Poisson Log Linear Mixture Models    #
# using EM to fit the model to the data.                  #
##=========================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("pmm-logLinear-EM.R")
library(R.utils)
library(edgeR)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)           # Set seed for reproducible results

##=============================================
# Generate the data and initialize variables  #
##=============================================
N         <- 1000               # Number of objects
sim       <- poisMixSim(N=N, libsize="A", separation="high")
X         <- sim$X              # Data points, Nxq matrix
conds     <- sim$conds          # Different biological conditions

K         <- 4                  # Number of clusters
epsilon   <- 1e-16              # Convergence paramater for EM
maxIter   <- 1000               # Maximum number of iterations for EM
eqProp    <- FALSE              # Equal mixing proportions pi
libSize   <- TRUE               # Take library size into consideration
libType   <- "TMM"              # Method to compute normalization factors
                                # Available: "TC", "MED", "DESEQ", "TMM"

# Wrap all the parameters in a list
theta     <- list(conds=conds, 
                  eqProp=eqProp,
                  libSize=libSize,
                  libType=libType)

##=====================================================
# Run PMM-LL-EM, explicitly giving initial parameters #
##=====================================================
fit.pmmLL <- pmm.LL.EM(X=X,
                       K=K,
                       theta=theta,
                       epsilon=epsilon,
                       maxIter=maxIter,
                       isDebug=TRUE)
