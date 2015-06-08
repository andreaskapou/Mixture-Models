##=======================================================================
# sumDerBin2ndPRL - Sum of the derivatives wrt to the parameters of     #
# the Binomial distributed Probit Regression Likelihood over all the    #
# dataset D, weighted by the corresponding responsibilites of each      #
# cluster. This function will be used in an Optimization procedure for  #
# the M-step of the EM algorithm                                        #
#                                                                       #
# Usage: res  <- sumBinomProbRegrLik(theta, D, post.resp)               #
#                                                                       #
# Input:                                                                #
#     theta   is a vector of the parameters of the function             #
#     D       is a list, where each entry contains a 3 by L dimensional #
#             matrix containing:                                        #
#         X.i is a vector (of length L) with the corresponding values   #
#               at each position for region i                           #
#         t.i is a vector (of length L) with the number of trials, in   #
#               the context of Binomal distributed data for region i    #
#         m.i is a vector (of length L) with the number of successes in #
#               the t trials for region i                               #
#     post.resp are the posterior probabilities, rersponsibility the    #
#               cluster k takes on data point x                         #
# Output:                                                               #
#     res     is a vector containing the total log likelihood of the    #
#               derivatives wrt to the parameters for all the regions i #
#                                                                       #
##=======================================================================
sumDerBin2ndPRL <- function(theta, D, post.resp){
  res <- c(0,0,0)
  for (i in 1:length(D)){
    res = res + bin2ndProbRegrLik(theta, D[[i]], mode=2) * post.resp[i]
  }
  return(res)
}