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


