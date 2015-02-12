logSumExp <- function(x) {
  ##=============================================
  # Function to normalize using the log sum exp #
  # trick for  numeric over/under flow          #
  ##=============================================
  # Computes log(sum(exp(x))
  offset <- max(x)
  return(log(sum(exp(x - offset))) + offset)
}