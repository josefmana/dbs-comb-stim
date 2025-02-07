#
# This is a script inspired by the DMC R package implementation (see osf.io/pbwx8)
# of statistical pipeline to analyse data from SSRT task under the assumption of
# two/three accumulator model with possibility of trigger & go failures
#


ssrt_simulate_participant <- function(
  
  n = 3,      # number of accumulators (first of which is a STOP accumulator)
  mu,         # a matrix with one row per accumulator, number of columns
  sigma,      # a vector of length one or the same format as mu
  tau,        # a vector of length one or the same format as mu
  TF = 0,     # probability of trigger failures
  GF = 0,     # probability of go failures
  SSD = Inf,  # vector of stop-signal delays
  Stair = 0.5 # adjustments via the staircase procedure

) {
  
}

