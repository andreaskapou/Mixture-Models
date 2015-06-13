#source("http://bioconductor.org/biocLite.R")
#biocLite("rowRanges")
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("binProbReg-EM.R")
library(R.utils)
library(ggplot2)
library(IRanges)
library(GenomicRanges)
library(BiSeq)
library(biomaRt)
library(matrixStats)
library(M3D)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(12345)                         # Set seed for reproducible results


#source("readENCODE.R")
#file      <- "GSM683770_hg19_wgEncodeHaibMethylRrbsH1hescHaibSitesRep1.bedRrbs.gz"
#methData  <- readENCODE(file)

data(rrbsDemo)
data(CpGsDemo)

overlaps  <- findOverlaps(CpGsDemo, rowData(rrbsDemo))

islands   <- unique(queryHits(overlaps))
CSites    <- rowData(rrbsDemo)
max       <- 172


X         <- list()
max <- 0
for (i in 1:length(islands)) {
  island <- islands[i]
  methIndices <- subjectHits(overlaps[queryHits(overlaps)==island])
  region <- CSites[methIndices]
  
  if (length(methIndices) > max){
    max <- length(methIndices)
    k <- i
  }
  meth  <- as.matrix(methReads(rrbsDemo)[methIndices,1])
  total <- as.matrix(totalReads(rrbsDemo)[methIndices,1])
  #values(region) <- data.frame(total=total, meth=meth)
  
  X[[i]] <- matrix(0, nrow=3, ncol=length(methIndices))
  #X[[i]][1,] <- seq(methIndices)/max
  X[[i]][1,] <- (methIndices-min(methIndices))/(max(methIndices)-min(methIndices))
  X[[i]][2,] <- total
  X[[i]][3,] <- meth
}

N       <- length(islands)
K       <- 5
deg     <- 4                              # Polynomial degree
epsilon <- 1e-4                           # Convergence paramater for EM
maxIter <- 100                            # Maximum number of iterations for EM

# Optional parameter initialization.
pi.c    <- rep(1/K, K)                    # Initialize mixing proportions for each cluster
theta   <- matrix(0, nrow=deg+1, ncol=K)  # Parameters of the 2nd order polynomial
for (k in 1:K){
  theta[,k] <- rep(0+k/100, deg+1)
}

params <- list(pi.c=pi.c, theta=theta)  # Wrap all the parameters in a list

##===========================================================
# Run BIN.PROB.REG-EM, explicitly giving initial parameters #
##===========================================================
fit.binProbReg <- binProbReg.EM(X=X, 
                                K=K, 
                                params=params, 
                                epsilon=epsilon, 
                                maxIter=22,
                                isDebug=TRUE)


##======================================================================
# Plot the results showing the K different functions that were learned #
##======================================================================
xs <- seq(0,1,len=2000) # create some values
# Add extra space to right of plot area; change clipping to figure
#par(xpd=T, mar=par()$mar+c(0,0,0,4))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

t <- 29
l <- 108
z <- 110

plot(X[[t]][1,], X[[t]][3,]/X[[t]][2,], col="darkgreen", pch=24, xlim=c(0,1), 
     ylim=c(0,1), xlab="region x", ylab="methylation level")
lines(x=xs, y=probPolynomFun(theta=fit.binProbReg$theta[,1], xs), 
      col="darkgreen", lwd=2)
points(X[[l]][1,], X[[l]][3,]/X[[l]][2,], pch=23, col="red")
lines(x=xs, y=probPolynomFun(theta=fit.binProbReg$theta[,3], xs), 
      col="red", lwd=2)
points(X[[z]][1,], X[[z]][3,]/X[[z]][2,], pch=8, col="darkblue")
lines(x=xs, y=probPolynomFun(theta=fit.binProbReg$theta[,5], xs), 
      col="darkblue", lwd=2)
