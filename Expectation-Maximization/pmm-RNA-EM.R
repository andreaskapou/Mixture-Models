##=========================================================
# Script for running Poisson Mixture Models using EM to   #
# fit the model to the data.                              #
##=========================================================

##=============================================================
# Set the working directory and load any required libraries   #
##=============================================================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
#source("pmm-RNA-EM.R")
library(R.utils)
library(edgeR)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)           # Set seed for reproducible results


N <- 200
K <- 4                    # Number of clusters
maxIter <- 1000
epsilon <- 1e-16
simulation  <- poisMixSim(N=N, libsize="A", separation="high")
X           <- simulation$X
conds       <- simulation$conds

lib.size    <- TRUE
lib.type    <- "TMM"

## Run the PMM model for g = 3
## "TC" library size estimate, EM algorithm

## Grouping columns of X in order of condition (all replicates put together)
o.ycols <- order(conds)
X       <- X[,o.ycols]
conds   <- conds[o.ycols]
conds.names <- unique(conds)    # Get unique condition names
d <- length(unique(conds))      # Total number of conditions
r <- as.vector(table(conds))    # Number of replicates in each condition
if(length(rownames(X)) == 0){   # If matrix X has no row names
  rn <- 1:nrow(X)
}
if(length(rownames(X)) > 0){ 
  rn <- rownames(X)
}

X <- as.matrix(X, nrow = NROW(X), ncol = NCOL(X))
rownames(X) <- rn               # Assign names to each row of X
N <- NROW(X)                    # Number of objects
q <- NCOL(X)                    # Number of variables
w <- rowSums(X)                 # Overall expression for each object

s <- rep(NA, q)                 # Normalized library size for each variable
if (lib.size == FALSE) {
  s <- rep(1, q)
}else{
  if (lib.type == "TC"){
    s <- colSums(X) / sum(X)
  }else if (lib.type == "MED"){
    s <- apply(X, MARGIN=2, FUN=median) / sum(apply(X, MARGIN=2, FUN=median))
  }else if (lib.type == "DESeq"){
    ## Code from DESeq function 'estimateSizeFactorsForMatrix'
    loggeomeans <- rowMeans(log(X))
    s <- apply(X, MARGIN=2, FUN=function(x) 
      exp(median((log(x)-loggeomeans)[is.finite(loggeomeans)])))
    s <- s / sum(s)
  }else if (lib.type == "TMM"){
    f <- calcNormFactors(X, method = "TMM")
    s <- colSums(X)*f / sum(colSums(X)*f)
  }
}

#s <- simulation$s.true
s.dot <- rep(NA, d)           # Sum of s on each replicate l
for (j in 1:d){
  s.dot[j] <- sum( s[which(conds == unique(conds)[j])] )
}

# Initialize parameters using 'kmeans'
cl    <- kmeans(X/w, K, nstart=25)        # Use Kmeans with random starts
C.n   <- cl$cluster                       # Get the mixture components
pi.c  <- as.vector(table(C.n)/NROW(X))    # Mixing proportions

C.mat <- matrix(0, nrow=N, ncol=K)        # Create matrix of cluster assignments
for(i in 1:N){
  C.mat[i, C.n[i]] <- 1
}

lambdas <- matrix(0, nrow=d, ncol=K)      # Calculate lambdas
denom   <- colSums(C.mat * w)
for(j in 1:d) {
  denom.norm  <- s.dot[j] * denom
  X.dot.j.dot <- rowSums(as.matrix(X[,which(conds == (unique(conds))[j])]))
  num         <- colSums( C.mat * matrix(rep(X.dot.j.dot, K), ncol=K))
  lambdas[j,] <- num / denom.norm
}

post.resp <- matrix(0, nrow=N, ncol=K)        # Hold responsibilities
pdf.w     <- matrix(0, nrow=N, ncol=K)        # Hold weighted PDFs
all.NLL   <- vector(mode="numeric")           # Hold NLL for all EM iterations
NLL       <- 1e+40                            # Initialize Negative Log Likelihood

mean.mat  <- vector("list", K)
w.mat     <- matrix(rep(w, times=q), nrow=N, ncol=q)
s.mat     <- matrix(rep(s, each=N) , nrow=N, ncol=q)

for (t in 1:maxIter){
  prevNLL   <- NLL                            # Store NLL to check for convergence
  
  for (k in 1:K){
    lambda.mat    <- matrix(rep(rep(lambdas[,k], times=r), each=N), nrow=N, ncol=q)
    mean.mat[[k]] <- w.mat * s.mat * lambda.mat
  }
  
  ##===================
  #       E-Step      #
  ##===================
  for (k in 1:K){
    pdf.w[,k] <- log(pi.c[k]) + rowSums(dpois(X, lambda=mean.mat[[k]], log=TRUE))
  }
  # Calculate probabilities using the logSumExp trick for numerical stability
  Z           <- apply(pdf.w, 1, logSumExp)
  post.resp   <- pdf.w - Z
  post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
  NLL <- -sum(Z)                            # Evaluate the NLL
  
  ##===================
  #       M-Step      #
  ##===================
  N.k     <- colSums(post.resp)             # Sum of responsibilities for each cluster
  pi.c    <- N.k/N                          # Update mixing proportions for each cluster
  
  denom   <- colSums(post.resp * w)
  for (j in 1:d){
    denom.norm  <- s.dot[j] * denom
    X.dot.j.dot <- rowSums(as.matrix(X[,which(conds == (unique(conds))[j])]))
    num         <- colSums(post.resp * matrix(rep(X.dot.j.dot, K), ncol=K))
    lambdas[j,] <- num / denom.norm
  }
  
  cat("NLL:", NLL, "\n")
  
  NLL.Diff  <- prevNLL - NLL                # Compute NLL difference after ith iteration
  if (NLL.Diff < 0){
    stop("Negative log likelihood increases - Something is wrong!")
  }
  all.NLL   <- c(all.NLL, NLL)              # Keep all NLL in a vector  
  if (NLL.Diff < epsilon){                  # Check for convergence.
    break
  }
}
