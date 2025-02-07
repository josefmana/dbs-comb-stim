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
    #"brms",      # for rexgaussian()
    #"cmdstanr",  # for model fitting
    "tidyverse", # for data wrangling
    "janitor",   # for 
    "bayesplot", # for MCMC-specific plotting
    "ggh4x"      # for special facet_wrapping
    
  )
)

# Load all in-house functions:
tar_source()

# Use multiple cores for model fitting:
#options( mc.cores = parallel::detectCores() )

# List the targets:
list(
  
  # CAUSAL ASSUMPTIONS ----
  tar_target(
    name    = DAG, # directed acyclic graph representing causal assumptions
    command = make_dag(plot = T, save = T)
  ),
  
  # FAKE DATA ----
  
  
  # MODEL(S) ----
  
  
  # REAL DATA IMPORT ----
  tar_target(
    name    = original_counterbalancing_path, # path to counterbalancing file
    command = here("_raw","rand","counterkey.csv"),
    format  = "file"
  ),
  tar_target(
    name    = original_counterbalancing_file, # read the file for counterbalancing
    command = read_delim(file = original_counterbalancing_path, delim = ";", escape_double = F, trim_ws = T)
  ),
  tar_target(
    name    = data_paths, # paths to single patients' data
    command = data_paths(folder = "ssrt")
  ),
  tar_target(
    name    = counterbalancing_data, # counterbalancing information
    command = extract_counterbalancing(.file = original_counterbalancing_file)
  ),
  tar_target(
    name    = demographic_file, # demographic data file
    command = here("_raw","demo","redcap_data.csv"),
    format  = "file"
  ),
  tar_target(
    name    = demographic_data, # extract demographic characteristics of included patients
    command = get_data(.paths = data_paths, .cb = counterbalancing_data, .demofile = demographic_file, type = "demography")
  ),
  tar_target(
    name    = ssrt_raw_data, # extract SSRT data (raw)
    command = get_data(.paths = data_paths, .cb = counterbalancing_data, .demofile = demographic_file, type = "raw")
  ),
  tar_target(
    name    = ssrt_labelled_data, # extract SSRT data (labelled)
    command = get_data(.paths = data_paths, .cb = counterbalancing_data, .demofile = demographic_file, type = "labelled")
  )
  
)
