
##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(R.utils)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory



fpr <- function(theta, n.i) {
  a   <- theta[1]
  b   <- theta[2]
  c   <- theta[3]
  m   <- n.i$m
  t   <- n.i$t
  x   <- m / t
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  res <- -dbinom(x=m, size=t, prob=Phi, log=TRUE)
  res <- sum(res)
  return(res)
}

gpr <- function(theta, n.i){
  a   <- theta[1]
  b   <- theta[2]
  c   <- theta[3]
  m   <- n.i$m
  t   <- n.i$t
  x   <- m / t
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  N   <- dnorm(g)
  
  da  <- sum(-N*x^2 * (m-t*Phi)/(Phi*(1-Phi)))
  db  <- sum(-N*x * (m-t*Phi)/(Phi*(1-Phi)))
  dc  <- sum(-N * (m-t*Phi)/(Phi*(1-Phi)))
  
  c(da, db, dc)
}


##====================
# Generate the data  #
##====================
epsilon <- 0.00001    # Convergence paramater
N       <- 200        # Total number of objects to create
K       <- 3          # Number of clusters
r1      <- rbinom(n=N, size=40, prob=.8)  # Total number of trials
r       <- matrix(r1, ncol=1)
X       <- gen.binomial(K=3, pi.c=c(.4,.3,.3), p=c(0.2,0.7, 0.5), r=r)

pdf.w       <- matrix(, N, K) # Hold the PDF of each point on each cluster k


for (k in 1:K){ # Calculate the PDF of each cluster for each data point
  pdf.w[,k] <- pi.c[k] * dbinom(X, size=r, prob=p[k])
}
post.resp   <- pdf.w / rowSums(pdf.w) # Get responsibilites by normalizarion

for (i in 1:1000){
  ## E-step
  for (k in 1:K) # Calculate the PDF of each cluster for each data point
    pdf.w[,k] <- log(pi.c[k]) + fpr(theta, n.i)
  post.resp   <-  pdf.w - apply(pdf.w, 1, logSumExp) # Normalize the log probability
  post.resp   <- apply(post.resp, 2, exp) # Eponentiate to get actual probabilities
  
  ## M-step
  for (k in 1:K){
    # Update mixing proportions
    pi.c[k]   <- sum(post.resp[, k] / sum(post.resp))
    
    p.optim <- optim(par=c(.6,.8,.2), fn=fpr, gr=gpr, n.i, method="CG")
}




t <- c(10, 11, 10, 11, 10)
m <- c(5, 7, 6, 3, 5)
x <- c(.7, .1, .2, .3, .8)
n.i <- list(m=m, t=t)
ss <- optim(par=c(.6,.8,.2), fn=fpr, gr=gpr, n.i, method="CG")