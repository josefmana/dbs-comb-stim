# This script is used to simulate SSRT data from the assumed data-generating process
# and fit Stan models with the same structure on the data

rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # set-up multiple cores

skip <- T # Skip model fitting if the model was already saved?

# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )
# install_cmdstan(version = "2.34.1")
# set_cmdstan_path()
#
# Important note: I was able to run the models successfully using cmdstanr version 2.34.1 on three different machines
# (one MacStudio and two MacBooks Pro, all with M1 or M2 processors), however, newer versions of cmdstanr stopped the chains
# due to the numerical integral not converging (using the same Stan code) and chains finishing unexpectedly!

library(here)
library(tidyverse)
#library(bayesplot)
library(cmdstanr)
library(ggh4x)

color_scheme_set("viridisA")
theme_set( theme_bw() )

if ( !dir.exists("sims") ) dir.create("sims") # folder for simulation results

source( here("scripts","fake_data_simulation.R") ) # read data generating function
source( here("scripts","utils.R") ) # utility functions

# In the next sections I generate data, fit models, do some posterior predictive checks, and parameter recovery checks.
# In the data used, there are 24 synthetic participants with realistic sample sizes (based on the real experiment) to get
# and idea of models' recovery properties and split them to chunks of eight (initial real sample size)
#
# Software used:
# R version 4.3.3 (2024-02-29)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.6.1)
#
# Since the code runs for quite some time, it is possible to read results generated on my laptop instead of running
# the analysis on yours (and the fitting section marked below by 'can be skipped' may be skipped):

if( file.exists( here("sims","ExGaussian_base.rds") ) ) fit0 <- readRDS( file = here("sims","ExGaussian_base.rds") )

# prepare the model
mod <- cmdstan_model( here("mods","ExGaussian_individual.stan") ) 

# function for manual initial values setting
ifun <- function() list(
  
  Int_mu_go_0 = runif(1,-2,0),
  Int_sigma_go_0 = runif(1,-2,0),
  Int_lambda_go_0 = runif(1,-2,0),
  Int_mu_stop_0 = runif(1,-2,0),
  Int_sigma_stop_0 = runif(1,-2,0),
  Int_lambda_stop_0 = runif(1,-2,0)
  
)


# MANY PARTICIPANTS, REALISTIC DATA ----

# generate data
# use default priors and default experiment length for 24 subjects
# setting priors of population-level variance to -Inf which translates to 0 on the log scale
# to ignore the hierarchical structure for now
d0 <- d0 <- ssrt_data_sim(

  tau_go = c(-Inf,0), tau_stop = c(-Inf,0),
  zeta_go = c(-Inf,0), zeta_stop = c(-Inf,0),
  epsilon_go = c(-Inf,0), epsilon_stop = c(-Inf,0),
  seeds = list(pars = 101:118, data = c(GO = 15, STOP = 68) ),
  N = 24

)

# some sanity checks
summary(d0$data$rt) # summaries
sum(d0$data$rt < 0, na.rm = T) / nrow(d0$data) # negative response times
t( sapply( unique(d0$data$id), function(i) table( subset(d0$data, signal == 1 & id == i)$response ) ) ) # success rates in STOP-SIGNAL trials (should be approximately 50/50 conditional on subject)
sanity_plot(d0$data) # GO response times should be generally slower than STOP-RESPOND response times

# loop through d0 data and fit the same individual model to each participant
# extract number of participants
k <- length( unique(d0$data$id) )

# skip of fit
if (skip == F) {
  
  # fit it
  fit0 <- lapply(
    
    X = set_names(x = 1:k),
    FUN = function(i) {
      
      current <- i # save current participant name for tracking
      
      # show a plot with current data being fitted highlighted 
      plot(
        
        sanity_plot(d0$data) + geom_rect(
          
          data = d0$data %>% mutate( subject = factor(paste0("Subject #",id), levels = unique( paste0("Subject #",id) ), ordered = T) ) %>% filter(id == current),
          aes(colour = subject, fill = "red"),
          xmin = -Inf, xmax = Inf,
          ymin = -Inf, ymax = Inf,
          alpha = 0.002,
          linewidth = 2.5
          
        )
        
      )
      
      # prepare data
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
      return( mod$sample(data = dlist, chains = 4, save_warmup = T, init = ifun, seed = 87542) )
      
    }
  )
  
  # save the results for immediate use if needed/wanted
  saveRDS( object = fit0, file = here("sims", "ExGaussian_base.rds") )
  
}


# plot trace plots, each for ca 3 s
for (i in 1:k) {

  print( paste0("participant #",i) )
  print( mcmc_trace( fit0[[i]]$draws(inc_warmup = T), n_warmup = 1e3 ) )
  Sys.sleep(3)

}


## ---- RECOVERY CHECKS ----

# prepare a data.frame with all parameter values (true and estimated)
pars <- left_join(
  
  # models' posteriors
  x = lapply(
    
    1:k, function(i) # loop through participants
      
      fit0[[i]]$draws(format = "data.frame") %>%
      select( starts_with("Int"), .chain, .iteration ) %>%
      mutate(ID = i, .before = 1)
    
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    rename_with( ~ sub( "Int_", "", sub("_0", "", .x) ) ) %>%
    mutate( # add estimated mean and SD from the parameters
      mean_go = exp(mu_go), sd_go = sqrt( exp(sigma_go)^2 + exp(lambda_go)^2 ),
      mean_stop = exp(mu_stop), sd_stop = sqrt( exp(sigma_stop)^2 + exp(lambda_stop)^2 )
    ) %>%
    pivot_longer(cols = contains("_"), values_to = "value", names_to = c("parameter", "type"), names_sep = "_"),
  
  # data-generating values
  y = with( d0$parameters, lapply(
    
    1:k, function(i)
      
      sapply(
        
        c("mu","sigma","lambda"), function(j) c(
          
          go = get(j)[type == "global intercept" & racer == "go"] + get(j)[type == "subject-level variability" & racer == "go"] * get(j)[type == "varying effect" & racer == "go" & id == i],
          stop = get(j)[type == "global intercept" & racer == "stop"] + get(j)[type == "subject-level variability" & racer == "stop"] * get(j)[type == "varying effect" & racer == "stop" & id == i]
          
        )
      ) %>%
      
      as.data.frame() %>%
      mutate( mean = exp(mu), sd = sqrt( exp(sigma)^2 + exp(lambda)^2) ) %>%
      rownames_to_column("type") %>%
      mutate(ID = i, .before = 1)
    
    )
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    pivot_longer( cols = -any_of( c("ID", "type") ), names_to = "parameter", values_to = "true_value" ),
  
  # glue estimated and true parameters to a single file
  by = c("ID","type","parameter")
  
) %>% left_join(
  
  # add observed mean and SDs for the go trials
  y = d0$data %>%
    filter(signal == 0) %>%
    group_by(id) %>%
    summarise( mean_go = mean(rt), sd_go = sd(rt) ) %>%
    ungroup() %>%
    pivot_longer(cols = contains("_"), values_to = "observed", names_to = c("parameter", "type"), names_sep = "_") %>%
    rename("ID" = "id"),
  
  # glue to the rest
  by = c("ID","type","parameter")
  
)

# check recovery of the parameters across participants
reco_hist( data = subset(pars, ID %in% 1:8), quants = c("mean", "sd") )
reco_hist( data = subset(pars, ID %in% 9:16), quants = c("mean", "sd") )
reco_hist( data = subset(pars, ID %in% 17:24), quants = c("mean", "sd") )


## ---- POSTERIOR PREDICTION ----

# extract a set of posterior predictions
ppred <- lapply(
  
  1:k, # one for each participant
  function(i) {
    
    # extract posterior draws for subject i
    df <-
      subset(pars, ID == i) %>%
      mutate( name = paste0(parameter,"_",type) ) %>%
      select(.chain, .iteration, name, value) %>%
      pivot_wider(names_from = name, values_from = value)
    
    sapply(
      
      1:nrow(df),
      function(j) {
        
        print( paste0("Participant #",i,", sample #",j) )
        
        return( with(
            
            df,
            ssrt_data_sim(
              alpha_go = c(mu_go[j], 0),
              alpha_stop = c(mu_stop[j], 0),
              beta_go = c(sigma_go[j], 0),
              beta_stop = c(sigma_stop[j], 0),
              gamma_go = c(lambda_go[j], 0),
              gamma_stop = c(lambda_stop[j], 0),
              tau_go = c(0, 0),
              tau_stop = c(0, 0),
              zeta_go = c(0, 0),
              zeta_stop = c(0, 0),
              epsilon_go = c(0, 0),
              epsilon_stop = c(0, 0),
              N = 1,
              df = subset(d0$data, id == i) %>% mutate(id = 1, rt = NA) %>% select(-ends_with("FT"), -winner)
            )$data$rt
          
          ) )
      }
      
    ) %>% t()

  }
)

# print posterior predictions (density overlaps)
# in the signal-respond trials, it makes little sense to "align" NAs between observed data and posterior predictions
ppc_density(d0$data, preds = ppred, cols = c("red4","blue4","lightpink2","skyblue2"), ncols = 4, ndrws = 50)

