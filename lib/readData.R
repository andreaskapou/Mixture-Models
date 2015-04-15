##=============================================
# Main script that will be called to read and #
# generate data                               #
##============================================= 

##=============================================
# Generate data:  mixture of 1D Gaussians     #
##=============================================
gen.gaussian <- function(N=500, K=2, pi.c=c(0.6,0.4), mus=c(10,20), stds=c(2,3)){
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Generate N random samples from a 1D-Gaussian, with the corresponding weights
  samples    <- rnorm(n=N,mean=mus[components],sd=stds[components])
  
  #plot(density(samples),main="Density Estimate of the Mixture Model")
  return (samples)
}

##=============================================
# Generate data:  mixture of MV Gaussians     #
##=============================================
gen.MV.gaussian <- function(N=500, K=2, pi.c=c(0.6,0.4), mus=matrix(c(1.5,4,4,1.5),2,2), S=list()){
  if (missing(S)){
    S[[1]] <- matrix(c(1,0,0,1),2,2)
    S[[2]] <- matrix(c(1,.7,.7,1),2,2)
  }
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Rows of the mean matrix are clusters and columns are the dimensionality of the data
  D         <- NCOL(mus)
  # Generate N random samples from a MV-Gaussian, with the corresponding weights
  samples   <- matrix(0, nrow=N, ncol=D)
  for (i in 1:N)
    samples[i,]    <- rmvnorm(1,mean=mus[components[i],], sigma=S[[components[i]]])
  return (samples)
}

##=============================================
# Generate data: mixture of 1D Poissons       #
##=============================================
gen.poisson <- function(N=500, K=2, pi.c=c(0.6,0.4), lambdas=c(8,22)){
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Generate N random samples from a Gaussian, with the corresponding weights
  samples    <- rpois(n=N,lambda=lambdas[components])
  
  return (samples)
}

##============================================
# Generate data: mixture of 1D Binomials     #
##============================================
gen.binomial <- function(N=500, K=2, pi.c=c(0.6,0.4), p=c(0.3,0.8), r){
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Generate N random samples from a Gaussian, with the corresponding weights
  samples    <- rbinom(n=N, size=r, prob=p[components])
  return (samples)
}

##=====================================
# Read data from Old Faithful Dataset #
##=====================================
read.faithful1D <- function(dim=2){
  X <-  faithful[,dim]
}
