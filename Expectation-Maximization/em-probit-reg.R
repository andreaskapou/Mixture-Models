
##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(1)


ff <- function(theta, X){
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  g <- pnorm(g)
  return(g)
}


l.theta <- function(theta, D, post.resp){
  res <- 0
  for (i in 1:length(D)){
    res = res + fpr(theta, D[[i]]) * post.resp[i]
  }
  return(res)
}

dl.theta <- function(theta, D, post.resp){
  res <- c(0,0,0)
  for (i in 1:length(D)){
    res = res + gpr(theta, D[[i]]) * post.resp[i]
  }
  return(res)
}

fpr <- function(theta, D){
  X   <- D[1,]
  t   <- D[2,]
  m   <- D[3,]
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  Phi <- pnorm(g)
  res <- sum(-dbinom(x=m, size=t, prob=Phi, log=TRUE))
  return(res)
}

gpr <- function(theta, D){
  X   <- D[1,]
  t   <- D[2,]
  m   <- D[3,]
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  Phi <- pnorm(g)
  N   <- dnorm(g)
  
  da  <- sum(-N*X^2 * (m-t*Phi)/(Phi*(1-Phi)))
  db  <- sum(-N*X * (m-t*Phi)/(Phi*(1-Phi)))
  dc  <- sum(-N * (m-t*Phi)/(Phi*(1-Phi)))
  
  c(da, db, dc)
}

##====================
# Generate the data  #
##====================
epsilon <- 1e-04      # Convergence paramater
N       <- 300        # Total number of objects to create
K       <- 3          # Number of clusters
pi.c    <- c(0.4, 0.3, 0.3) # Mixing proportions

theta   <- matrix(0, nrow=3, ncol=K)
for (k in 1:K){
  theta[,k] <- c(0+k/100, 0+k/100, 0+k/100)
}

theta[,1] <- c(-0.024, -4, -1.01)
theta[,2] <- c(0.024, 4, -1.01)
theta[,3] <- c(-4, 0, 0.4)

X       <- gen.meth.data2(N=N, pi.c=pi.c) # Generate methylation profiles data

pdf.w       <- matrix(, N, K) # Hold the PDF of each point on each cluster k
post.resp   <- matrix(, N, K)
log.likel   <- 0

for (t in 1:10){
  ## E-step
  for (k in 1:K){ # Calculate the PDF of each cluster for each data point
    for (i in 1:N){
      pdf.w[i,k] <- log(pi.c[k]) + fpr(theta=theta[,k], X[[i]])
    }
  }
  
  post.resp   <- pdf.w - apply(pdf.w, 1, logSumExp) # Normalize the log probability
  post.resp   <- apply(post.resp, 2, exp) # Exponentiate to get actual probabilities
  
  ## M-step
  for (k in 1:K){
    # Update mixing proportions
    pi.c[k]   <- sum(post.resp[,k] / sum(post.resp))
    
    opt.theta <- optim(par=theta[,k], fn=l.theta, gr=dl.theta, X, post.resp[,k], method="CG", 
                                                                    control = list(maxit = 10))
    theta[,k] <- opt.theta$par
  }
}

#theta[,1] <- c(-0.024, -4, -1.01)
#theta[,1] <- c(-4, 0, 0.4)
#theta[,1] <- c(0.024, 4, -1.01)
xs <- seq(-1,1,len=100) # create some values
plot(x=xs, y=ff(theta=theta[,1], xs), type="l", xlab="x", ylab="", 
     main=expression(a*x^2 + b*x + c))
