# This script is used to simulate SSRT data from the assumed.

rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # set-up multiple cores

# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )

library(here)
library(tidyverse)
library(brms) # for rexgaussian()
library(bayesplot)
library(cmdstanr)

color_scheme_set("viridisA")
theme_set( theme_bw() )


# DATA SIMULATION FUNCTION ----

ssrt_data_sim <- function( alpha_go = c(-0.4,0.2), # global intercept of the go racer mu parameter
                           alpha_stop = c(-1.0,0.2), # global intercept of the stop racer mu parameter
                           beta_go = c(-2.0,0.2), # global intercept of the go racer sigma parameter
                           beta_stop = c(-2.0,0.2), # global intercept of the stop racer sigma parameter
                           gamma_go = c(-2.0,0.2), # global intercept of the go racer lambda parameter
                           gamma_stop = c(-2.0,0.2), # global intercept of the stop racer lambda parameter
                           tau_go = c(-2.0,0.2), # subject-level standard deviation of the go racer mu parameter
                           tau_stop = c(-2.0,0.2), # subject-level standard deviation of the stop racer mu parameter
                           zeta_go = c(-2.0,0.2), # subject-level standard deviation of the go racer sigma parameter
                           zeta_stop = c(-2.0,0.2), # subject-level standard deviation of the stop racer sigma parameter
                           epsilon_go = c(-2.0,0.2), # subject-level standard deviation of the go racer lambda parameter
                           epsilon_stop = c(-2.0,0.2), # subject-level standard deviation of the stop racer lambda parameter
                           N = 8, # number of subjects
                           K = c(216,72) # number of go/stop-signal trials
                        ) {
  
  
  ## SAMPLE EXGAUSSIAN PARAMETERS ----
  
  # sample global intercepts
  alphaGO <- rnorm( 1, alpha_go[1], alpha_go[2] )
  alphaSTOP <- rnorm( 1, alpha_stop[1], alpha_stop[2] )
  betaGO <- rnorm( 1, beta_go[1], beta_go[2] )
  betaSTOP <- rnorm( 1, beta_stop[1], beta_stop[2] )
  gammaGO <- rnorm( 1, gamma_go[1], gamma_go[2] )
  gammaSTOP <- rnorm( 1, gamma_stop[1], gamma_stop[2] )
  
  # sample subject-level standard deviations
  tauGO <- exp( rnorm( 1, tau_go[1], tau_go[2] ) )
  tauSTOP <- exp( rnorm( 1, tau_stop[1], tau_stop[2] ) )
  zetaGO <- exp( rnorm( 1, zeta_go[1], zeta_go[2] ) )
  zetaSTOP <- exp( rnorm( 1, zeta_stop[1], zeta_stop[2] ) )
  epsilonGO <- exp( rnorm( 1, epsilon_go[1], epsilon_go[2] ) )
  epsilonSTOP <- exp( rnorm( 1, epsilon_stop[1], epsilon_stop[2] ) )
  
  # sample standardised subject-level effects
  # if only one participant is generated, ignore the subject level part by setting these to zero
  for( i in c("xGO","yGO","zGO","xSTOP","ySTOP","zSTOP") ) if (N == 1) assign( i, 0 ) else assign( i, rnorm(N) )
  
  # calculate GO and STOP ExGaussian parameters
  muGO = exp(alphaGO + xGO * tauGO)
  muSTOP = exp(alphaSTOP + xSTOP * tauSTOP)
  sigmaGO = exp(betaGO + yGO * zetaGO)
  sigmaSTOP = exp(betaSTOP + ySTOP * zetaSTOP)
  lambdaGO = exp(gammaGO + zGO * epsilonGO)
  lambdaSTOP = exp(gammaSTOP + zSTOP * epsilonSTOP)
  
  
  ## CONDUCT A VIRTUAL EXPERIMENT ----
  
  # prepare the order of experimental coniditions
  # each subject receives a column with K[2] trial numbers with a stop signal
  signal_trials <- sapply( 1:N, function(i) sample(1:sum(K), K[2], replace = F, prob = NULL) )
  
  # prepare the output data matrix
  out <- lapply(

    1:N,
    function(i)
      matrix( data = NA, nrow = sum(K), ncol = 6, dimnames = list( trial = NULL, variable = c("trial","id","signal","ssd","response","rt") ) )
  
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
  
  # loop through all the rows of the output matrix
  for ( i in 1:nrow(out) ) {
    
    ### GO TRIALS ----

    if( out[ i, "signal"] == 0 ) {
      
      out[ i, "response" ] <- 1 # assume correct response
      out[ i, "rt" ] <- rexgaussian(1, muGO[out[i,"id"]], sigmaGO[out[i,"id"]], lambdaGO[out[i,"id"]] ) # sample response time
      
    ### STOP TRIALS ----
    
    } else if( out[ i, "signal"] == 1 ) {
      
      # set-up initial SSD if it is subject's first stop-signal trial
      if( out[i, "trial"] == min( signal_trials[ , out[i,"id"]] ) ) SSD <- .3
      
      # write down SSD for this trial
      out[ i, "ssd" ] <- SSD
      
      # sample finishing times of GO and STOP racers
      goFT <- rexgaussian(1, muGO[out[i,"id"]], sigmaGO[out[i,"id"]], lambdaGO[out[i,"id"]] )
      stopFT <- SSD + rexgaussian(1, muSTOP[out[i,"id"]], sigmaSTOP[out[i,"id"]], lambdaSTOP[out[i,"id"]] )
      
      # fill-in the rest depending on the winner
      # if GO racer wins
      if (goFT < stopFT) {
        
        out[ i, "response" ] <- 1 # incorrect response is recorded
        out[ i, "rt"] <- goFT # with GO racer finishing time as the response time
        SSD <- SSD - .05 # make it easier during the next trial
        
      # else if STOP racer wins
      } else if (stopFT < goFT) {
        
        out[ i, "response" ] <- 0 # correct non-response with no response time is recorded
        SSD <- SSD + .05 # make it harder during the next trial
      }
    }

  }
  
  ## PREPARE THE OUTCOMES ----

  # ExGaussian racers' parameters
  pars <- data.frame(
    
    id = c( rep(NA,4), rep(1:N, 2) ),
    type = c( rep( c("global intercept","subject-level variability"), 2 ), rep("varying effect", 2*N) ),
    racer = c( rep("go", 2), rep("stop", 2), rep("go", N), rep("stop", N) ),
    mu = c(alphaGO, tauGO, alphaSTOP, tauSTOP, xGO, xSTOP),
    sigma = c(betaGO, zetaGO, betaSTOP, zetaSTOP, yGO, ySTOP),
    lambda = c(gammaGO, epsilonGO, gammaSTOP, epsilonSTOP, zGO, zSTOP)

  )
  
  # re-format data for later tidyverse shinaningans
  dats <- as.data.frame(out)
  
  # return list of parameters and data
  return( list( parameters = pars, data = dats ) )
  
}


# GENERATE DATA ----

# use default priors and default experiment length for 8 subjects
d0 <- ssrt_data_sim(N = 8)

## SOME SANITY CHECKS ----

# check for negative response times
summary( d0$data$rt )
sum(d0$data$rt < 0, na.rm = T) / nrow(d0$data)

# check success rates in STOP-SIGNAL trials
# should be approximately 50/50 conditional on subject
t( sapply( unique(d0$data$id), function(i) table( subset(d0$data, signal == 1 & id == i)$response ) ) )

# GO response times should be generally slower than STOP-RESPOND response times
# under the assumption of the horse race model
d0$data %>%
  
  # some formatting shinanigans
  filter( complete.cases(rt) ) %>%
  mutate( `Response type: ` = ifelse( signal == "0", "GO", "SIGNAL-RESPOND" ) ) %>%
  mutate( subject = paste0("Subject #",id) ) %>%
  
  # plot it
  ggplot() +
  aes(x = rt, fill = `Response type: `) +
  geom_density(alpha = .5) +
  labs(x = "Response time (s)", y = "Density") +
  facet_wrap(~subject, ncol = 2) +
  theme(legend.position = "bottom")


# FIT MODELS ----

## INDIVIUDAL EXGAUSSIAN MODEL ----

# prepare the model
mod_indi <- cmdstan_model( here("mods","ExGaussian_individual.stan") )

# function for manual initial values setting
ifun <- function() list(

  Int_mu_go_0 = runif(1,-2,0),
  Int_sigma_go_0 = runif(1,-2,0),
  Int_lambda_go_0 = runif(1,-2,0),
  Int_mu_stop_0 = runif(1,-2,0),
  Int_sigma_stop_0 = runif(1,-2,0),
  Int_lambda_stop_0 = runif(1,-2,0)

)


### single participant realistic data set ----

# this chunk uses simulate data of a single participant
# in an experiment as big as our real data to check
# convergence in realistically sized data
d1 <- ssrt_data_sim(N = 1)

# prepare separate files for "go" and "signal-respond", and "successful inhibition" trials
Dgo <- d1$data %>% filter( signal == 0 )
Dsr <- d1$data %>% filter( signal == 1 & !is.na(rt) )
Dna <- d1$data %>% filter( signal == 1 & is.na(rt) )

# prepare input data
dlist1 <- list(
  
  # data
  Y_0_go = Dgo$rt, N_0_go = nrow(Dgo),
  Y_0_sr = Dsr$rt, N_0_sr = nrow(Dsr), SSD_0_sr = Dsr$ssd,
  N_0_na = nrow(Dna), SSD_0_na = Dna$ssd,
  
  # priors
  mu_go_0_p = c(-0.4,0.2), sigma_go_0_p = c(-2.0,0.2), lambda_go_0_p = c(-2.0,0.2),
  mu_stop_0_p = c(-1.0,0.2), sigma_stop_0_p = c(-2.0,0.2), lambda_stop_0_p = c(-2.0,0.2)
  
)

# model fitting
fit1 <- mod_indi$sample( data = dlist1, chains = 4, save_warmup = T, init = ifun )

# results checking
#mcmc_trace( fit1$draws(inc_warmup = T), n_warmup = 1e3 ) # trace plots
fit1$summary() # parameter estimates
d1$parameters # true parameters


### single participant bulky data set ----

# this chunk uses simulate data of a single participant
# in an experiment much bigger than our real data to check
# convergence in 'large' numbers
#d2 <- ssrt_data_sim(K = c(3e3,1e3), N = 1)

# prepare separate files for "go" and "signal-respond", and "successful inhibition" trials
#Dgo <- d2$data %>% filter( signal == 0 )
#Dsr <- d2$data %>% filter( signal == 1 & !is.na(rt) )
#Dna <- d2$data %>% filter( signal == 1 & is.na(rt) )

# prepare input data
#dlist2 <- list(
#  
#  # data
#  Y_0_go = Dgo$rt, N_0_go = nrow(Dgo),
#  Y_0_sr = Dsr$rt, N_0_sr = nrow(Dsr), SSD_0_sr = Dsr$ssd,
#  N_0_na = nrow(Dna), SSD_0_na = Dna$ssd,
#  
#  # priors
#  mu_go_0_p = c(-0.4,0.2), sigma_go_0_p = c(-2.0,0.2), lambda_go_0_p = c(-2.0,0.2),
#  mu_stop_0_p = c(-1.0,0.2), sigma_stop_0_p = c(-2.0,0.2), lambda_stop_0_p = c(-2.0,0.2)
#  
#)

# model fitting
#fit2 <- mod_indi$sample( data = dlist2, chains = 4, save_warmup = T, init = ifun )

# results checking
#mcmc_trace( fit2$draws(inc_warmup = T), n_warmup = 1e3 ) # trace plots
#fit2$summary() # parameter estimates
#d2$parameters # true parameters


### single model, many participants ----

# loop through d0 data and fit the same individual model to each participant
# extract number of participants
k <- length( unique(d0$data$id) )

# prepare a folder for the fits
fit_many_indi <- list()

# fit it
for ( i in 1:k ) {
  
  # data
  dGO <- subset(d0$data, id == i) %>% filter( signal == 0 )
  dSR <- subset(d0$data, id == i) %>% filter( signal == 1 & !is.na(rt) )
  dNA <- subset(d0$data, id == i) %>% filter( signal == 1 & is.na(rt) )
  
  # input file
  dlist <- list(
    
    # data
    Y_0_go = dGO$rt, N_0_go = nrow(dGO),
    Y_0_sr = dSR$rt, N_0_sr = nrow(dSR), SSD_0_sr = dSR$ssd,
    N_0_na = nrow(dNA), SSD_0_na = dNA$ssd,
    
    # priors
    mu_go_0_p = c(-0.4,0.2), sigma_go_0_p = c(-2.0,0.2), lambda_go_0_p = c(-2.0,0.2),
    mu_stop_0_p = c(-1.0,0.2), sigma_stop_0_p = c(-2.0,0.2), lambda_stop_0_p = c(-2.0,0.2)
    
  )
  
  # fitting proper
  fit_many_indi[[i]] <- mod_indi$sample( data = dlist, chains = 4, save_warmup = T, init = ifun )
  
}

# print median parameter estimates
lapply(
  
  1:k,
  function(i)
    t( fit_many_indi[[i]]$summary()[ , c("variable","median") ] %>% column_to_rownames("variable") )
  
) %>%
  
  do.call( rbind, . ) %>%
  `rownames<-`(1:k)

# print the true values
with(
  
  d0$parameters,
  lapply(
    
    setNames( c("mu","sigma","lambda"), c("mu","sigma","lambda") ),
    function(i)
      
      sapply(
        1:k,
        function(j)
          c( go = get(i)[type == "global intercept" & racer == "go"] + get(i)[type == "subject-level variability" & racer == "go"] * get(i)[type == "varying effect" & racer == "go" & id == j],
             stop = get(i)[type == "global intercept" & racer == "stop"] + get(i)[type == "subject-level variability" & racer == "stop"] * get(i)[type == "varying effect" & racer == "stop" & id == j]
             )
        
      )
  )
)


## HIERARCHICAL EXGAUSSIAN MODEL ----

# prepare the model
mod_hier <- cmdstan_model( here("mods","ExGaussian_hierarchical.stan") )

# extract number of participants
k <- length( unique(d0$data$id) )

# function for manual initial values setting
ifun_hier <- function() list(
  
  alpha_go_0 = runif(1,-2,0), beta_go_0 = runif(1,-2,0), gamma_go_0 = runif(1,-2,0),
  alpha_stop_0 = runif(1,-2,0), beta_stop_0 = runif(1,-2,0), gamma_stop_0 = runif(1,-2,0),
  
  tau_go_0 = runif(1,-2,0), zeta_go_0 = runif(1,-2,0), epsilon_go_0 = runif(1,-2,0),
  tau_stop_0 = runif(1,-2,0), zeta_stop_0 = runif(1,-2,0), epsilon_stop_0 = runif(1,-2,0),
  
  x_go_0 = rep(0,k), y_go_0 = rep(0,k), z_go_0 = rep(0,k),
  x_stop_0 = rep(0,k), y_stop_0 = rep(0,k), z_stop_0 = rep(0,k)
  
)


### realistically sized data set ----

# let us use the d0 file simulates before
# prepare separate files for "go" and "signal-respond", and "successful inhibition" trials
Dgo <- d0$data %>% filter( signal == 0 )
Dsr <- d0$data %>% filter( signal == 1 & !is.na(rt) )
Dna <- d0$data %>% filter( signal == 1 & is.na(rt) )

# prepare input data
dlist3 <- list(
  
  # data
  N_0_go = nrow(Dgo), Y_0_go = Dgo$rt, J_0_go = Dgo$id,
  N_0_sr = nrow(Dsr), Y_0_sr = Dsr$rt, SSD_0_sr = Dsr$ssd, J_0_sr = Dsr$id,
  N_0_na = nrow(Dna), SSD_0_na = Dna$ssd, J_0_na = Dna$id,
  K = k,
  
  # group-level priors
  alpha_go_0_p = c(-0.4,0.2), beta_go_0_p = c(-2.0,0.2), gamma_go_0_p = c(-2.0,0.2),
  alpha_stop_0_p = c(-1.0,0.2), beta_stop_0_p = c(-2.0,0.2), gamma_stop_0_p = c(-2.0,0.2),
  
  # participant-level priors
  tau_go_0_p = c(-2.0,0.2), zeta_go_0_p = c(-2.0,0.2), epsilon_go_0_p = c(-2.0,0.2),
  tau_stop_0_p = c(-2.0,0.2), zeta_stop_0_p = c(-2.0,0.2), epsilon_stop_0_p = c(-2.0,0.2)
  
)

# model fitting
fit3 <- mod_hier$sample( data = dlist3, chains = 4, save_warmup = T, init = ifun_hier )

# results checking
#mcmc_trace( fit3$draws(inc_warmup = T), n_warmup = 1e3 ) # trace plots
fit3$summary() # parameter estimates
d0$parameters # true parameters

