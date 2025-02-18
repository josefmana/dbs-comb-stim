#
# This script defines a function that simulates SSRT data based on two racers model.
#

ssrt_data_simulate <- function(
  
  alpha_go = c(-0.4,0.2),     # global intercept of the go racer mu parameter
  alpha_stop = c(-1.0,0.2),   # global intercept of the stop racer mu parameter
  beta_go = c(-2.0,0.2),      # global intercept of the go racer sigma parameter
  beta_stop = c(-2.0,0.2),    # global intercept of the stop racer sigma parameter
  gamma_go = c(-2.0,0.2),     # global intercept of the go racer lambda parameter
  gamma_stop = c(-2.0,0.2),   # global intercept of the stop racer lambda parameter
  tau_go = c(-2.0,0.2),       # subject-level standard deviation of the go racer mu parameter
  tau_stop = c(-2.0,0.2),     # subject-level standard deviation of the stop racer mu parameter
  zeta_go = c(-2.0,0.2),      # subject-level standard deviation of the go racer sigma parameter
  zeta_stop = c(-2.0,0.2),    # subject-level standard deviation of the stop racer sigma parameter
  epsilon_go = c(-2.0,0.2),   # subject-level standard deviation of the go racer lambda parameter
  epsilon_stop = c(-2.0,0.2), # subject-level standard deviation of the stop racer lambda parameter
  multilevel = F,             # use multilevel specification?
  N = 8,                      # number of subjects
  K = c(216,72),              # number of go/stop-signal trials
  S = NULL,                   # seeds, either a vector of length 18 or NULL for parameters
  df = NULL                   # a fixed data set for posterior predictive checks
  
) {
  
  ## SAMPLE EXGAUSSIAN PARAMETERS ----
  
  # prepare seeds
  if ( is.null(S) ) S <- sample(x = 1:1e9, size = 18, replace = F)
  
  # make a matrix of it
  S <- matrix(

    data = S,
    ncol = 2,
    byrow = T,
    dimnames = list(
      x = c("alpha","beta","gamma","tau","zeta","epsilon","x","y","z"),
      y = c("GO","STOP")
    )

  )
  
  # How many intercepts to draw from the commont distribution?
  # Single intercept per type for a multilevel specification (where variability in parameters comes from subject-level variance)
  # N different intercepts for a single-level specification (where varibility in parameters comes from prior only)
  n <- ifelse( multilevel == T, 1, N)
  
  # sample global intercepts
  set.seed(S["alpha","GO"]); alphaGO <- rnorm( n, alpha_go[1], alpha_go[2] )
  set.seed(S["alpha","STOP"]); alphaSTOP <- rnorm( n, alpha_stop[1], alpha_stop[2] )
  set.seed(S["beta","GO"]); betaGO <- rnorm( n, beta_go[1], beta_go[2] )
  set.seed(S["beta","STOP"]); betaSTOP <- rnorm( n, beta_stop[1], beta_stop[2] )
  set.seed(S["gamma","GO"]); gammaGO <- rnorm( n, gamma_go[1], gamma_go[2] )
  set.seed(S["gamma","STOP"]); gammaSTOP <- rnorm( n, gamma_stop[1], gamma_stop[2] )
  
  
  # sample subject-level standard deviations
  set.seed(S["tau","GO"]); tauGO <- ifelse(multilevel == T, exp( rnorm( 1, tau_go[1], tau_go[2] ) ), 0)
  set.seed(S["tau","STOP"]); tauSTOP <- ifelse(multilevel == T, exp( rnorm( 1, tau_stop[1], tau_stop[2] ) ), 0)
  set.seed(S["zeta","GO"]); zetaGO <- ifelse(multilevel == T, exp( rnorm( 1, zeta_go[1], zeta_go[2] ) ), 0)
  set.seed(S["zeta","STOP"]); zetaSTOP <- ifelse(multilevel == T, exp( rnorm( 1, zeta_stop[1], zeta_stop[2] ) ), 0)
  set.seed(S["epsilon","GO"]); epsilonGO <- ifelse(multilevel == T, exp( rnorm( 1, epsilon_go[1], epsilon_go[2] ) ), 0)
  set.seed(S["epsilon","STOP"]); epsilonSTOP <- ifelse(multilevel == T, exp( rnorm( 1, epsilon_stop[1], epsilon_stop[2] ) ), 0)
  
  # sample standardised subject-level effects
  # if only one participant is generated, ignore the subject level part by setting these to zero
  for( i in c("xGO","yGO","zGO","xSTOP","ySTOP","zSTOP") ) {
    
    if (N == 1) assign( i, 0 )
    else {
      set.seed(S[substr(i, 1, 1), substr(i, 2, nchar(i))])
      assign( i, rnorm(N) )
    }
    
  }
  
  # calculate GO and STOP ExGaussian parameters
  muGO = exp(alphaGO + xGO * tauGO)
  muSTOP = exp(alphaSTOP + xSTOP * tauSTOP)
  sigmaGO = exp(betaGO + yGO * zetaGO)
  sigmaSTOP = exp(betaSTOP + ySTOP * zetaSTOP)
  lambdaGO = exp(gammaGO + zGO * epsilonGO)
  lambdaSTOP = exp(gammaSTOP + zSTOP * epsilonSTOP)
  
  
  ## CONDUCT A VIRTUAL EXPERIMENT ----
  
  ### ---- PREPARE DATASET ----
  
  # if there is no template for posterior predictive checks
  # prepare the data from ground
  if ( is.null(df) ) {
    
    # prepare the order of experimental conditions
    # each subject receives a column with K[2] trial numbers with a stop signal
    signal_trials <- sapply( 1:N, function(i) sample(1:sum(K), K[2], replace = F, prob = NULL) )
    
    # prepare the output data matrix
    out <- lapply(
      
      1:N,
      function(i)
        matrix(
          data = NA,
          nrow = sum(K),
          ncol = 3,
          dimnames = list(
            trial = NULL,
            variable = c("trial","id","signal")
          )
        )
      
    )
    
    # pre-allocate
    for( i in 1:N ) {
      
      out[[i]][ , "trial"] <- 1:sum(K) # trial numbers
      out[[i]][ , "id"] <- i # IDs
      out[[i]][ signal_trials[ ,i], "signal"] <- 1 # stop-signal trials
      out[[i]][ (1:sum(K))[ !(1:sum(K)) %in% signal_trials[ ,i] ], "signal"] <- 0 # go trials
      
    }
    
    # pull into a single file
    out <- do.call(rbind, out)
    
  } else {
    
    out <- df[ , c("trial", "id", "signal") ]
    signal_trials <- sapply( 1:N, function(i) which( subset(df, id == i)$signal == 1) )
    
  }
  
  
  ### ---- FILL-IN RESPONSES ----
  
  out <- cbind(
    
    out,
    sapply(
      
      1:nrow(out), function(i) {
        
        #### ---- GO TRIALS ----
        
        if( out[i, "signal"] == 0 ) {
          
          response <- 1 # assume correct response
          rt <- rexgaussian(1, muGO[out[i,"id"]], sigmaGO[out[i,"id"]], lambdaGO[out[i,"id"]] ) # sample response time
          SSD <- NA
          
          #### ---- STOP TRIALS ----
        } else if( out[i, "signal"] == 1 ) {
          
          # set-up initial SSD if it is subject's first stop-signal trial
          if( out[i, "trial"] == min( signal_trials[ , out[i,"id"]] ) ) SSD <- .3 else SSD <- SSD_next
          
          # sample finishing times of GO and STOP racers
          goFT <- rexgaussian(1, muGO[out[i,"id"]], sigmaGO[out[i,"id"]], lambdaGO[out[i,"id"]] )
          stopFT <- SSD + rexgaussian(1, muSTOP[out[i,"id"]], sigmaSTOP[out[i,"id"]], lambdaSTOP[out[i,"id"]] )
          
          
          # fill-in the rest depending on the winner
          # if GO racer wins
          if (goFT < stopFT) {
            
            response <- 1 # incorrect response is recorded
            rt <- goFT # with GO racer finishing time as the response time
            SSD_next <<- SSD - .05 # make it easier during the next trial
            
            # else if STOP racer wins
          } else if (stopFT < goFT) {
            
            response <- 0 # correct non-response with no response time is recorded
            rt <- NA
            SSD_next <<- SSD + .05 # make it harder during the next trial
          }
        }
        
        return( c(response = response, rt = rt, ssd = SSD) )
      }
    ) %>% t()
    
  )


  ## PREPARE THE OUTPUT ----
  
  if (multilevel == T) {
    
    pars <- data.frame(
      
      id = c( rep(NA,4), rep(1:N, 2) ),
      type = c( rep( c("global intercept","subject-level variability"), 2 ), rep("varying effect", 2*N) ),
      racer = c( rep("go", 2), rep("stop", 2), rep("go", N), rep("stop", N) ),
      mu = c(alphaGO, tauGO, alphaSTOP, tauSTOP, xGO, xSTOP),
      sigma = c(betaGO, zetaGO, betaSTOP, zetaSTOP, yGO, ySTOP),
      lambda = c(gammaGO, epsilonGO, gammaSTOP, epsilonSTOP, zGO, zSTOP)
      
    )
    
  } else pars <- data.frame(
    
    id = c( rep(NA,4), rep(1:N, 2) ),
    type = c( rep( c("global intercept","subject-level variability"), 2 ), rep("varying effect", 2*N) ),
    racer = c( rep("go", 2), rep("stop", 2), rep("go", N), rep("stop", N) ),
    mu = c( rep(c(0,1),2), alphaGO, alphaSTOP ),
    sigma = c( rep(c(0,1),2), betaGO, betaSTOP ),
    lambda = c( rep(c(0,1),2), gammaGO, gammaSTOP )
    
  )
  
  # re-format data for later tidyverse shinanigans
  dats <- as.data.frame(out)
  
  # return list of parameters and data
  return( list( seeds = S, parameters = pars, data = dats ) )
  
}

