# Dirichlet Process sampling (Stick-breaking construction, with unbounded number of sticks)
# 
# Input
# - alpha: dispersion/ concentration parameter
# - H0: a baseline measure, a function H0() which gives samples 
#
# Output
# - a discrete measure sampled from H0 (pi_1, pi_2, ...), with its associated values (theta_1, theta_2, ...)
#
# Stick-breaking
# - pi_k = beta * prod_{i<k}(pi_i)
# - beta ~ Beta(1, alpha)
#
# Alternative construction: fix the number of sticks to a parameter k.
rdirprocess = function(alpha, H0)
{
	require(MCMCpack)
	
	# Mixture proportions (need sum_i pi_i = 1). In the non-parametric formulation these are unbounded in number 
	pi = NULL 

 	# Values of the proportions sampled from H0, for instance the means of Gaussian clusters. Same size as pi
 	values = NULL
	
 	# this will be sum_i pi_i
 	culsums = 0
	
	# Stick-breaking: we go as far as the sum of the mixture proprotions pi integrates to 1.
	while(sum(pi) < 1)
	{
		# ~ H0: draw a new sample from the baseline distribution
		values = c(values, H0()) 
			
		# the new proportion pi_k = prod_{j<n} pi_j * Beta(1, alpha)
		pi = c(pi, prod(1-pi) * rbeta(1,1,alpha)) 
	}
	
	# renormalize pi becase we might have slightly gone over 1
	pi = pi/sum(pi)
	return(list(pi=pi, values=values))	
}

# Random deviate from empirical disitribution (here to sample from the disitribution drawn from DP(alpha,H0))
rempd = function(n, pi ,values)
{
	samples = NULL
	for(i in 1:n)
	{
		r = runif(1)
		
	idx = which(cumsum(pi) <= r)
	if(length(idx) == 0) idx = 1
	else idx = max(idx)
	
		samples  = c(samples, values[idx])
	}
	
	return(samples)
}

# higher-order function for a Gaussian sampling -- H0(0,1) is a sampler for N(0,1)
H0 = function(mu, sigma){ return(function() rnorm(1, mu, sigma)) }

# We make some plots for some values of the dispersion parameter alpha
dev.new(height = 18)
par(mfrow = c(7,2))
for(alpha in c(0.1, 1, 10, 10, 100, 1000, 10000))
{
	# plot ~ N(0,1) rescaled in [0,1]
	x   <- seq(-5,+5,length=100)
	y   <- dnorm(x,mean=0, sd=1)
	y = y/max(y)
	plot(x,y, type="l", lwd=2)

	# ~ DP(alpha H0) sample, plot with probabilities rescaled in [0,1]
	x = rdirprocess(alpha, H0(0,1))
	lines(x = x$values, y = x$pi/max(x$pi), col = 'red', pch=16, type = 'h')
	
	h = hist(x$values, plot = F)
	h$density = h$counts/max(h$counts)
	plot(h,freq=FALSE, add = TRUE, border = 'darkgreen')

	title(paste('alpha = ', alpha))
	
	# 10K sample from the empirical distribution
	synthdata = rempd(10000, x$pi, x$values)
	h = hist(synthdata, plot = F)
	h$density = h$counts/max(h$counts)
	plot(h, freq=FALSE)
	lines(density(synthdata), col = 'brown', lwd = 2)
}
dev.copy2pdf(file = 'DP-sample-example.pdf')

