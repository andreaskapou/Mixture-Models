
##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory

ff <- function(theta, X){
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  g <- pnorm(g)
  return(g)
}


fpr <- function(theta, D){
  X   <- D[1]
  t   <- D[2]
  m   <- D[3]
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  Phi <- pnorm(g)
  res <- sum(-dbinom(x=m, size=t, prob=Phi, log=TRUE))
  return(res)
}

gpr <- function(theta, D){
  X   <- D[1]
  t   <- D[2]
  m   <- D[3]
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  Phi <- pnorm(g)
  N   <- dnorm(g)
  
  da  <- sum(-N*X^2 * (m-t*Phi)/(Phi*(1-Phi)))
  db  <- sum(-N*X * (m-t*Phi)/(Phi*(1-Phi)))
  dc  <- sum(-N * (m-t*Phi)/(Phi*(1-Phi)))
  
  c(da, db, dc)
}


#t <- c(10, 11, 11, 11, 10, 10)
#m <- c(8, 7, 7, 3, 2, 2)
#X <- c(-50, -40, -30, 20, 30, 50)/(50/2)
#ss <- optim(par=c(0,0,0), fn=fpr, gr=gpr, X, t, m, method="CG")


# plot the function
#xs <- seq(-1,1,len=100) # create some values
#plot(x=xs, y=ff(theta=ss$par, xs), type="l", xlab="x", ylab="", 
#     main=expression(a*x^2 + b*x + c))

##====================
# Generate the data  #
##====================
epsilon <- 0.00001    # Convergence paramater
N       <- 300        # Total number of objects to create
K       <- 3          # Number of clusters
pi.c    <- c(0.45, 0.35, 0.2) # Mixing proportions

theta   <- matrix(0, nrow=K, ncol=3)
for (i in 1:K){
  theta[i,] <- c(0+i/100, 0+i/100, 0+i/100)
}

X       <- gen.meth.data(N=N, pi.c=pi.c) # Generate methylation profiles data

pdf.w       <- matrix(, N, K) # Hold the PDF of each point on each cluster k
post.resp   <- matrix(, N, K)
log.likel   <- 0

for (t in 1:1000){
  ## E-step
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    for (i in 1:N){
      pdf.w[i,k] <- log(pi.c[k]) + fpr(theta=theta[k,], X[[i]])
    }
  }
  post.resp   <-  pdf.w - apply(pdf.w, 1, logSumExp) # Normalize the log probability
  post.resp   <- apply(post.resp, 2, exp) # Exponentiate to get actual probabilities
  
  ## M-step
  for (k in 1:K){
    # Update mixing proportions
    pi.c[k]   <- sum(post.resp[, k] / sum(post.resp))
    
    theta.optim <- optim(par=theta[k,], fn=fpr, gr=gpr, X[[i]], method="CG")
  }
}

