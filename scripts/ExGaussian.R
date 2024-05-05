# This script is used to diagnose ExGaussian model of SSRT data implemented in Stan.

rm( list = ls() )
options( mc.cores = parallel::detectCores() )

# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )

library(here)
library(tidyverse)
#library(brms)
library(bayesplot)
library(cmdstanr)

# read data
d0 <- read.csv( here("_data","ssrt_lab.csv"), sep = "," )

# pivot it wider
d1 <-
  d0 %>%
  filter( block != 0 ) %>%
  select( id, block, cond, signal, rt1, trueSOA ) %>%
  rename( "rt" = "rt1", "ssd" = "trueSOA" )

# separate files for "go" and "stop" trials
Dgo <- d1 %>% filter( signal == "nosignal" )
Dsr <- d1 %>% filter( signal == "signal" & !is.na(rt) )
Dna <- d1 %>% filter( signal == "signal" & is.na(rt) )


# PRIOR PREDICTION ----

# prepare a function computing and summarising mean and sd of ExGaussian based on priors
priorpred <- function( mu = c(700,50), sigma = c(0,2), tau = c(1.7,1.3), sd = c(0,100), rep = 1e5 ) {
  
  m <- rnorm( rep, mu[1], mu[2] ) + rnorm( rep, 0, abs( rnorm( rep, sd[1], sd[2] ) ) ) + exp( rnorm( rep, tau[1], tau[2] ) )
  sd <- sqrt( exp( rnorm( rep, sigma[1], sigma[2] ) )^2 + exp( rnorm( rep, tau[1], tau[2] ) )^2 )
  
  return( summary( cbind(M = m, SD = sd) ) )
  
}

# MODEL FITTING ----

# drop rows with missing data in Go trials
GOcon <- subset( Dgo, complete.cases(rt) & cond == "ctrl" ) %>% mutate( id = as.integer( as.factor(id) ) )
GOexp <- subset( Dgo, complete.cases(rt) & cond == "exp" ) %>% mutate( id = as.integer( as.factor(id) ) )

# prepare the model
mod <- cmdstan_model( here("mods","ExGaussian_hier.stan") )

# prepare input
dlist <- list(
  
  # data
  Y_0 = GOcon$rt, N_0 = nrow(GOcon), J_0 = GOcon$id, Z_0_1 = rep( 1, nrow(GOcon) ),
  Y_1 = GOexp$rt, N_1 = nrow(GOexp), J_1 = GOexp$id, Z_1_1 = rep( 1, nrow(GOexp) ),
  K = length( unique(GOcon$id) ), M = 1,
  
  # priors
  mu_0p = c(700,50), sigma_0p = c(1,2), beta_0p = c(1,2),
  mu_1p = c(700,50), sigma_1p = c(1,2), beta_1p = c(1,2),
  tau_p = c(0,200)
  
)

# fit the model
fit <- mod$sample( data = dlist, chains = 4 )

