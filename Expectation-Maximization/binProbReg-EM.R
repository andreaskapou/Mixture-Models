
##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(R.utils)
library(ggplot2)
sourceDirectory("../lib", modifiedOnly=FALSE) # Source the 'lib' directory
set.seed(4)

##====================
# Generate the data  #
##====================
epsilon     <- 1e-04            # Convergence paramater
N           <- 400              # Total number of objects to create
K           <- 3                # Number of clusters
pi.c        <- c(.45, .35, .2)  # Mixing proportions
theta       <- matrix(0, nrow=3, ncol=K) # Parameters of the 2nd order polynomial
for (k in 1:K){
  theta[,k] <- c(0+k/100, 0+k/100, 0+k/100)
}
X           <- gen.meth.data(N=N, pi.c=pi.c) # Generate methylation profiles

pdf.w       <- matrix(, N, K)   # Hold the PDF of each point on each cluster k
post.resp   <- matrix(, N, K)   # Hold responsibilities
logLik      <- 0

for (t in 1:100){               # Loop until convergence
  prevLogLik  <- logLik         # Store to check for convergence
  
  ##========
  # E-Step #
  ##========
  for (k in 1:K){ # Calculate the weighted PDF of each cluster for each data point
    for (i in 1:N){
      pdf.w[i,k] <- log(pi.c[k]) + binomProbRegrLik(theta=theta[,k], D=X[[i]], mode=1)
    }
  }
  # ln p(X|mu,S,p) = Sum_{n=1}^{N}(ln(Sum_{k=1}^{K}(p_k * Bin(x_n|Phi(g)))))
  logLik      <- sum(log(colSums(exp(pdf.w))))      # Evaluate the log likelihood
  
  post.resp   <- pdf.w - apply(pdf.w, 1, logSumExp) # Normalize the log probability
  post.resp   <- apply(post.resp, 2, exp)           # Exponentiate to get actual probabilities
  
  ##========
  # M-Step #
  ##========
  for (k in 1:K){
    N.k       <- sum(post.resp[,k])                 # Sum of responsibilities for cluster k
    pi.c[k]   <- N.k / N                            # Update mixing proportions for cluster k
    
    # Update the parameters a, b, c of the polynomial using Conjugate-Gradient method
    theta[,k] <- optim(par=theta[,k],
                       fn=sumBinomProbRegrLik,
                       gr=sumDerBinomProbRegrLik, X, post.resp[,k],
                       method="CG",
                       control = list(fnscale=-1, maxit = 5) )$par
  }
  if (abs(logLik - prevLogLik) < epsilon){          # Check for convergence.
    break
  }
  print(logLik)
}


##=================================================
# Simple function for second order polynomial     #
# transformed through through the probit          #
# function, and thus it is squashed to be in the  #
# (0,1) interval                                  #
##=================================================
ff <- function(theta, X){
  g   <- theta[1]*X^2 + theta[2]*X + theta[3]
  g   <- pnorm(g)
  return(g)
}
xs <- seq(-1,1,len=2000) # create some values

plot(x=xs, y=ff(theta=theta[,1], xs), type="l", xlab="x", ylab="", 
     main=expression(a*x^2 + b*x + c), col="darkgreen")
lines(x=xs, y=ff(theta=theta[,2], xs), col="red")
lines(x=xs, y=ff(theta=theta[,3], xs), col="blue")


par(xpd=T, mar=par()$mar+c(0,0,0,4))
# Add extra space to right of plot area; change clipping to figure
#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(X[[7]][1,], X[[7]][3,]/X[[7]][2,], col="darkgreen", pch=24, xlim=c(-1,1), ylim=c(0.05,0.92),
     xlab="region x", ylab="methylation level")
lines(x=xs, y=ff(theta=theta[,1], xs), col="darkgreen", lwd=2)
points(X[[245]][1,], X[[245]][3,]/X[[245]][2,], pch=23, col="red")
lines(x=xs, y=ff(theta=theta[,3], xs), col="red", lwd=2)
points(X[[390]][1,], X[[390]][3,]/X[[390]][2,], pch=8, col="darkblue")
lines(x=xs, y=ff(theta=theta[,2], xs), col="darkblue", lwd=2)

# Add legend to top right, outside plot region
legend("right", inset=c(-0.1,0), legend=c("1","2", "3"), pch=c(24,23,8), 
       col=c("darkgreen", "red", "darkblue"), title="Cluster")


#legend(locator(1),c("Cluster", "group B"), pch = c(1,2), lty = c(1,2))

dev.copy(png, width = 800, height = 600, filename=paste("../images/probitParams", format(Sys.time(), "%a%b%d%H%M"),".png", sep=""))
dev.off()

#save(theta, file=paste("../files/probitTheta", format(Sys.time(), "%a%b%d%H%M"),".RData", sep=""))
