#Finite Mixture Models using Gibbs Sampling

Code for modeling data using finite Dirichlet mixture models under a Bayesian approach.
It is implemented for Gaussian mixtures, Binomial mixtures and Poisson mixtures. 

1. For the mixing proportions a Dirichlet conjugate prior is used. 
2. For the Gaussian mixture a NormalGamma conjuagate prior is set for the mean and the precision
3. For the Binomial mixture a Beta conjugate prior is set for the success probabilities
4. For the Poisson mixture a Gamma conjugate prior is set for the mean/variance parameter

By using conjugate priors, the Gibbs sampling procedure can be used to update the full conditionals
on each MCMC iteration. 
