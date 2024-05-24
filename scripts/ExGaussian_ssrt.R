# This script is used to run ExGaussian model of SSRT data implemented in Stan.

rm( list = ls() )
options( mc.cores = parallel::detectCores() )

# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )

library(here)
library(tidyverse)
library(brms) # for rexgaussian()
library(bayesplot)
library(cmdstanr)

# maximum time allow was 1.5 s
rtmax <- 1.5

# read data
d0 <- read.csv( here("_data","ssrt_lab.csv"), sep = "," )

# pivot it wider
d1 <-
  d0 %>%
  mutate( cens = ifelse( correct == "missed", "right", "none" ), max = rtmax ) %>% # add censoring info and maximum possible RTs
  filter( block != 0 ) %>%
  select( id, block, cond, signal, rt1, trueSOA, cens, max ) %>%
  rename( "rt" = "rt1", "ssd" = "trueSOA" ) %>%
  mutate( across( c("rt","ssd"), ~ .x/1e3 ) ) # re-scale from ms to s

# separate files for "go" and "stop" trials
Dgo <- d1 %>% filter( signal == "nosignal" )
Dsr <- d1 %>% filter( signal == "signal" & !is.na(rt) )
Dna <- d1 %>% filter( signal == "signal" & is.na(rt) ) %>% mutate( rt = rtmax, cens = "right" )


# MODEL FITTING ----

# prepare data sets for model
GOcon <- subset( Dgo, complete.cases(rt) & cond == "ctrl" ) %>% mutate( id = as.integer( as.factor(id) ) )
GOexp <- subset( Dgo, complete.cases(rt) & cond == "exp" ) %>% mutate( id = as.integer( as.factor(id) ) ) 
SRcon <- subset( Dsr, complete.cases(rt) & cond == "ctrl" & (rt>ssd) ) %>% mutate( id = as.integer( as.factor(id) ) )
SRexp <- subset( Dsr, complete.cases(rt) & cond == "exp" & (rt>ssd) ) %>% mutate( id = as.integer( as.factor(id) ) )
NAcon <- subset( Dna, cond == "ctrl" ) %>% mutate( id = as.integer( as.factor(id) ) )
NAexp <- subset( Dna, cond == "exp" ) %>% mutate( id = as.integer( as.factor(id) ) )

# prepare the model
mod <- cmdstan_model( here("mods","ExGaussian_ssrt.stan") )

# prepare input
dlist <- list(
  
  # data
  Y_0_go = GOcon$rt, N_0_go = nrow(GOcon), J_0_go = GOcon$id, Z_0_go = rep( 1, nrow(GOcon) ),
  Y_1_go = GOexp$rt, N_1_go = nrow(GOexp), J_1_go = GOexp$id, Z_1_go = rep( 1, nrow(GOexp) ),
  Y_0_sr = SRcon$rt, N_0_sr = nrow(SRcon), J_0_sr = SRcon$id, Z_0_sr = rep( 1, nrow(SRcon) ), SSD_0_sr = SRcon$ssd,
  Y_1_sr = SRexp$rt, N_1_sr = nrow(SRexp), J_1_sr = SRexp$id, Z_1_sr = rep( 1, nrow(SRexp) ), SSD_1_sr = SRexp$ssd,
  #Y_0_na = NAcon$rt, N_0_na = nrow(NAcon), J_0_na = NAcon$id, Z_0_na = rep( 1, nrow(NAcon) ), SSD_0_na = NAcon$ssd,
  #Y_1_na = NAexp$rt, N_1_na = nrow(NAexp), J_1_na = NAexp$id, Z_1_na = rep( 1, nrow(NAexp) ), SSD_1_na = NAexp$ssd,
  
  K = length( unique(GOcon$id) ), M = 1,
  
  # priors
  mu_go_0_p = c(-.4,.2), sigma_go_0_p = c(-2,.2), beta_go_0_p = c(-2,.2),
  mu_go_1_p = c(-.4,.2), sigma_go_1_p = c(-2,.2), beta_go_1_p = c(-2,.2),
  mu_stop_0_p = c(-.4,.2), sigma_stop_0_p = c(-2,.2), beta_stop_0_p = c(-2,.2),
  mu_stop_1_p = c(-.4,.2), sigma_stop_1_p = c(-2,.2), beta_stop_1_p = c(-2,.2),
  tau_mu_go_p = c(0,.5), tau_sigma_go_p = c(0,.3), tau_beta_go_p = c(0,.3),
  tau_mu_stop_p = c(0,.5), tau_sigma_stop_p = c(0,.3), tau_beta_stop_p = c(0,.3)
  
)

# fit the model
fit <- mod$sample( data = dlist, chains = 4 )

# check trace plots
mcmc_trace( fit$draws() )


# POSTERIOR PREDICTION ----

# extract draws
post <- fit$draws(format = "matrix")

# predict separately results for control and experimental condition
d_seq <- with(

  dlist,
  list( con = data.frame( rt = Y_0, id = J_0, cond = "con" ), exp = data.frame( rt = Y_1, id = J_1, cond = "exp" ) )
  
)

# predict control condition data
con_pred <-
  
  sapply(
    
    1:nrow(d_seq$con),
    function(i) {
      
      print( paste0("row #",i," out of ",nrow(d_seq$con) ) )
      
      rexgaussian(
        nrow(post),
        mu = exp( post[ ,"InterceptMu_0"] + post[ ,paste0("r_mu[",d_seq$con$id[i],"]") ] ),
        sigma = exp( post[ ,"InterceptSigma_0"] + post[ ,paste0("r_sigma[",d_seq$con$id[i],"]") ] ),
        beta = exp( post[ ,"InterceptBeta_0"] + post[ ,paste0("r_beta[",d_seq$con$id[i],"]") ] )
      )

    }
  )

# predict experimental condition data
exp_pred <-
  
  sapply(
    
    1:nrow(d_seq$exp),
    function(i) {
      
      print( paste0("row #",i," out of ",nrow(d_seq$exp) ) )
      
      rexgaussian(
        nrow(post),
        mu = exp( post[ ,"InterceptMu_1"] + post[ ,paste0("r_1_mu[",d_seq$exp$id[i],"]") ] ),
        sigma = exp( post[ ,"InterceptSigma_1"] + post[ ,paste0("r_1_sigma[",d_seq$exp$id[i],"]") ] ),
        beta = exp( post[ ,"InterceptBeta_1"] + post[ ,paste0("r_1_beta[",d_seq$exp$id[i],"]") ] )
      )
      
    }
  )

# check posterior stats
d_pred <-
  do.call( rbind.data.frame, d_seq ) %>%
  mutate( conid = paste0(cond,"_",id) ) %>%
  mutate( conid = factor( conid, levels = c( paste0("con_",1:4), paste0("exp_",1:4), paste0("con_",5:8), paste0("exp_",5:8) ), ordered = T ) )

ppc_stat_grouped( y = d_pred$rt, yrep = cbind(con_pred,exp_pred), group = d_pred$conid, stat = "mean" )
ppc_dens_overlay_grouped( y = d_pred$rt, yrep = cbind(con_pred,exp_pred)[ sample(1:4000,100) , ], group = d_pred$conid )
