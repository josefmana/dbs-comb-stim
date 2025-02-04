#
# This is a script running targets pipeline of the project.
#

# Load packages required to define the pipeline:
library(targets)
#library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c(
    
    "here",      # for path listing
    "ggdag",     # for causal diagrams
    "ggraph",    # for graphs 
    "brms",      # for rexgaussian()
    "tidyverse", # for data wrangling
    "janitor",   # for 
    "cmdstanr",  # for model fitting
    "bayesplot", # for MCMC-specific plotting
    "ggh4x"      # for special facet_wrapping
    
  )
)

# Load all in-house functions:
tar_source()

# Use multiple cores for model fitting:
options( mc.cores = parallel::detectCores() )


# List the targets:
list(
  
  ## CAUSAL ASSUMPTIONS ----
  tar_target(
    name = DAG, # directed acyclic graph representing causal assumptions
    command = make_dag(plot = T, save = T)
  ),
  
  ## MODEL(S) ----
  tar_target(
    name = stan_model, # the Stan model to be used
    command = here("ExGaussian_individual.stan"),
    format = "file"
  ),
  tar_target(
    name = stat_model, # read the model
    command = cmdstan_model(stan_model)
  ),
  
  # FAKE DATA ----
  tar_target(
    name = fake_data, # simulate
    command = ssrt_data_simulate(N = 8, multilevel = F, S = 49:66)
  ),
  tar_target(
    name = sanity_check, # check for negative response times
    command = data_summary(fake_data$data, 5)
  ),
  tar_target(
    name = sanity_plot, # check observed data distributions
    command = response_times_plot(fake_data$data, ncols = 2)
  ),
  
    ### ---- model fitting ----
    tar_target(
      name = fit,
      command = fit_individually(fake_data$data, stat_model)
    ),
    tar_target(
      name = trace_plots,
      command = show_trace(fit)
    ),
  
    ### ---- recovery checks ----
    tar_target(
      name = quantities,
      command = extract_parameters(fit, fake_data$parameters, fake_data$data)
    ),
    tar_target(
      name = recovery_plot_means,
      command = recovery_plot( data = quantities, quants = c("mean", "sd") )
    ),
    tar_target(
      name = recovery_plot_parameters,
      command = recovery_plot( data = quantities, quants = c("mu", "sigma", "lambda"), tit = "Parameters estimates" )
    ),
  
    ### ---- posterior predictive check ----
    tar_target(
      name = posterior_predictions,
      command = compute_predictions(fit, quantities, fake_data$data)
    ),
    tar_target(
      name = ppc_density_plot,
      command = ppc_density(fake_data$data, posterior_predictions, cols = c("red4","blue4","lightpink2","skyblue2"), ncols = 2, ndrws = 50)
    ),
  
  ## COUNTERBALANCING ----
  tar_target(
    name = original_counterbalancing_path, # path to counterbalancing file
    command = here("_raw","rand","counterkey.csv"),
    format = "file"
  ),
  tar_target(
    name = original_counterbalancing_file, # read the file for counterbalancing
    command = read_delim(file = original_counterbalancing_path, delim = ";", escape_double = F, trim_ws = T)
  ),
  tar_target(
    name = retest_counterbalancing_options, # generate valid options of counterbalancing in re-test
    command = generate_counterbalancing(.file = original_counterbalancing_file)
  ),
  
  ## DATA IMPORT ----
  tar_target(
    name = data_paths, # paths to single patients' data
    command = data_paths(folder = "ssrt")
  ),
  tar_target(
    name = counterbalancing_data, # counterbalancing information
    command = extract_counterbalancing(.file = original_counterbalancing_file)
  ),
  tar_target(
    name = demographic_file, # demographic data file
    command = here("_raw","demo","redcap_data.csv"),
    format = "file"
  ),
  tar_target(
    name = demographic_data, # extract demographic characteristics of included patients
    command = get_data(.paths = data_paths, .cb = counterbalancing_data, .demofile = demographic_file, type = "demography")
  ),
  tar_target(
    name = ssrt_raw_data, # extract SSRT data (raw)
    command = get_data(.paths = data_paths, .cb = counterbalancing_data, .demofile = demographic_file, type = "raw")
  ),
  tar_target(
    name = ssrt_labelled_data, # extract SSRT data (labelled)
    command = get_data(.paths = data_paths, .cb = counterbalancing_data, .demofile = demographic_file, type = "labelled")
  )
  
)
