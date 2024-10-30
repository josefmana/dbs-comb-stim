# This is a script running targets pipeline of the project.

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set( packages = c(

  "here", # for path listing
  "ggdag", # for causal diagrams
  "ggraph", # for graphs 
  "brms", # for rexgaussian()
  "tidyverse", # for data wrangling
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
  
  # CAUSAL ASSUMPTIONS ----
  tar_target( DAG, make_dag(plot = F, save = T) ),
  
  # MODEL(S) ----
  tar_target( stan_model, "ExGaussian_individual.stan", format = "file" ), # the Stan model to be used
  tar_target( stat_model, cmdstan_model(stan_model) ), # read the model
  
  # FAKE DATA ----
  tar_target( fake_data, ssrt_data_simulate(N = 8, multilevel = F, S = 49:66) ), # simulate
  tar_target( sanity_check, data_summary(fake_data$data, 5) ), # check for negative response times
  tar_target( sanity_plot, response_times_plot(fake_data$data, ncols = 2) ), # check observed data distributions
  
    ## ---- model fitting ----
    tar_target( fit, fit_individually(fake_data$data, stat_model) ),
    tar_target( trace_plots, show_trace(fit) ),
  
    ## ---- recovery checks ----
    tar_target( quantities, extract_parameters(fit, fake_data$parameters, fake_data$data) ),
    tar_target( recovery_plot_means, recovery_plot( data = quantities, quants = c("mean", "sd") ) ),
    tar_target( recovery_plot_parameters, recovery_plot( data = quantities, quants = c("mu", "sigma", "lambda"), tit = "Parameters estimates" ) ),
  
    ## ---- posterior predictive check ----
    tar_target( posterior_predictions, compute_predictions(fit, quantities, fake_data$data) ),
    tar_target( ppc_density_plot, ppc_density(fake_data$data, posterior_predictions, cols = c("red4","blue4","lightpink2","skyblue2"), ncols = 2, ndrws = 50) )

)
