# This is a script running targets pipeline of the project.

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set( packages = c(

  "here", # for path listing
  "brms", # for rexgaussian()
  "tidyverse", # for data wranngling
  "cmdstanr", # for model fitting
  "bayesplot", # for MCMC-specific plotting
  "ggh4x" # for special facet_wrapping

) )

# Load all in-house functions:
tar_source()

# Use multiple cores for model fitting:
options( mc.cores = parallel::detectCores() )


# List the targets:
list(
  
  # Fake data simulation ----
  
  tar_target(stan_model, "ExGaussian_individual.stan", format = "file"), # the Stan model to be used
  tar_target(mod0, cmdstan_model(stan_model) ), # read the model
  
  # simulate data 
  tar_target(d0, ssrt_data_sim( N = 8, multilevel = F, S = 57:74) ),
  
  # run some sanity checks
  tar_target( san_check, fake_data_sums(d0$data, 5) ),
  tar_target( san_plot, sanity_plot(d0$data, ncols = 2) ),
  
  # fit the models
  tar_target( fit0, indi_fit(d0$data, mod0) ),
  
  # check convergence
  tar_target( trace0, show_trace(fit0) ),
  
  # recovery checks
  tar_target( pars0, get_pars(fit0, d0$parameters, d0$data) ),
  tar_target( recoplot_locscale, reco_hist( data = pars0, quants = c("mean", "sd") ) ),
  tar_target( recoplot_paramets, reco_hist( data = pars0, quants = c("mu", "sigma", "lambda"), tit = "Parameters estimates" ) ),
  
  # posterior predictive check
  tar_target( ppred, ppred_calc(fit0, pars0, d0$data) ),
  tar_target( ppc_dens, ppc_density(d0$data, ppred, cols = c("red4","blue4","lightpink2","skyblue2"), ncols = 2, ndrws = 50) )

)
