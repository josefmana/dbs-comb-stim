#
# This script is used to simulate SSRT data from the assumed data-generating process
# and fit Stan models with the same structure on the data
#
# recommended running this is a fresh R session or restarting current session
# install.packages( "cmdstanr", repos = c( "https://mc-stan.org/r-packages/", getOption("repos") ) )
# install_cmdstan(version = "2.34.1")
# set_cmdstan_path()
#
# Important note: I was able to run the models successfully using cmdstanr version 2.34.1 on three different machines
# (one MacStudio and two MacBooks Pro, all with M1 or M2 processors), however, newer versions of cmdstanr stopped the chains
# due to the numerical integral not converging (using the same Stan code) and chains finishing unexpectedly!
#
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
#

#
# ---- MANUAL INITIAL VALUES SETTINGS ----
ifun <- function() list(
  
  Int_mu_go = runif(1,-2,0),
  Int_sigma_go = runif(1,-2,0),
  Int_lambda_go = runif(1,-2,0),
  Int_mu_stop = runif(1,-2,0),
  Int_sigma_stop = runif(1,-2,0),
  Int_lambda_stop = runif(1,-2,0)
  
)

#
# ---- SANITY CHECKS ----
data_summary <- function(data, wait) with(
  
  data, {
    
    print( summary(rt) ) # summaries
    print( sum(rt < 0, na.rm = T) / length(rt) ) # negative response times
    print( t( sapply( unique(id), function(i) table( response[signal == 1 & id == i] ) ) ) ) # success rates in STOP-SIGNAL trials (should be approximately 50/50 conditional on subject)
    Sys.sleep(wait) # stop for a second so that we can check it.
    
  }
)

#
# ---- INDIVIDUAL MODELS FITTING ----
fit_individually <- function(data, model) {
  
  # extract number of participants
  k <- length( unique(data$id) )
  
  # fit it
  fit <- lapply(
    
    X = set_names(x = 1:k),
    FUN = function(i) {
      
      current <- i # save current participant name for tracking
      
      # FOR USE OUT OF THE TARGETS PIPELINE
      # show a plot with current data being fitted highlighted 
      #plot(
      #  
      #  sanity_plot(data) + geom_rect(
      #    
      #    data = data %>% mutate( subject = factor(paste0("Subject #",id), levels = unique( paste0("Subject #",id) ), ordered = T) ) %>% filter(id == current),
      #    aes(colour = subject, fill = "red"),
      #    xmin = -Inf, xmax = Inf,
      #    ymin = -Inf, ymax = Inf,
      #    alpha = 0.002,
      #    linewidth = 2.5
      #    
      #  )
      #  
      #)
      
      # prepare data
      dGO <- subset(data, id == i) %>% filter( signal == 0 )
      dSR <- subset(data, id == i) %>% filter( signal == 1 & !is.na(rt) )
      dNA <- subset(data, id == i) %>% filter( signal == 1 & is.na(rt) )
      
      # input file
      dlist <- list(
        
        # data
        Y_go = dGO$rt, N_go = nrow(dGO),
        Y_sr = dSR$rt, N_sr = nrow(dSR), SSD_sr = dSR$ssd,
        N_na = nrow(dNA), SSD_na = dNA$ssd,
        
        # priors
        prior_muGO = c(-0.4,0.2), prior_sigmaGO = c(-2.0,0.2), prior_lambdaGO = c(-2.0,0.2),
        prior_muSTOP = c(-1.0,0.2), prior_sigmaSTOP = c(-2.0,0.2), prior_lambdaSTOP = c(-2.0,0.2),
        
        # integrand & debugging specification
        integration_mode = 1, integrand_variant = 0, debug_integrand = 1
        
      )
      
      # prepare a folder for the outcome
      folder <- here( "_sims", paste0("indi_fake_subject_no_",i) )
      if( dir.exists(folder) ) unlink(folder, recursive = T)
      dir.create(folder)
      
      # fitting proper
      return( model$sample(
        
        data = dlist,
        chains = 4,
        save_warmup = T,
        init = ifun,
        seed = 87542,
        output_dir = folder
        
      ) )
      
    }
    
  )
  
  # return it
  return(fit)

}

#
# ---- TRACE PLOTS ----
show_trace <- function(fit) lapply(
  
  X = 1:length(fit),
  FUN = function(i) mcmc_trace( fit[[i]]$draws(inc_warmup = T), n_warmup = 1e3 )

)

#
# ---- EXTRACT MODEL PARAMETERS ----
extract_parameters <- function(fit, truth, data) left_join(
  
  # models' posteriors
  x = lapply(
    
    1:length(fit), function(i) # loop through participants
      
      fit[[i]]$draws(format = "data.frame") %>%
      select( starts_with("Int"), .chain, .iteration ) %>%
      mutate(ID = i, .before = 1)
    
  ) %>%
    
    do.call( rbind.data.frame, . ) %>%
    rename_with( ~ sub( "Int_", "", .x) ) %>%
    mutate( # add estimated mean and SD from the parameters
      mean_go = exp(mu_go), sd_go = sqrt( exp(sigma_go)^2 + exp(lambda_go)^2 ),
      mean_stop = exp(mu_stop), sd_stop = sqrt( exp(sigma_stop)^2 + exp(lambda_stop)^2 )
    ) %>%
    pivot_longer(cols = contains("_"), values_to = "value", names_to = c("parameter", "type"), names_sep = "_"),
  
  # data-generating values
  y = with( truth, lapply(
    
    1:length(fit), function(i)
      
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
    pivot_longer(cols = -any_of( c("ID", "type") ), names_to = "parameter", values_to = "true_value"),
  
  # glue estimated and true parameters to a single file
  by = c("ID","type","parameter")
  
) %>% left_join(
  
  # add observed mean and SDs for the go trials
  y = data %>%
    filter(signal == 0) %>%
    group_by(id) %>%
    summarise( mean_go = mean(rt), sd_go = sd(rt) ) %>%
    ungroup() %>%
    pivot_longer(cols = contains("_"), values_to = "observed", names_to = c("parameter", "type"), names_sep = "_") %>%
    rename("ID" = "id"),
  
  # glue to the rest
  by = c("ID","type","parameter")
  
)

#
# ---- POSTERIOR PREDICTION ----
compute_predictions <- function(fit, pars, data) lapply(
  
  1:length(fit), # one for each participant
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
          ssrt_data_simulate(
            alpha_go = c(mu_go[j], 0), alpha_stop = c(mu_stop[j], 0),
            beta_go = c(sigma_go[j], 0), beta_stop = c(sigma_stop[j], 0),
            gamma_go = c(lambda_go[j], 0), gamma_stop = c(lambda_stop[j], 0),
            tau_go = c(-Inf, 0), tau_stop = c(-Inf, 0),
            zeta_go = c(-Inf, 0), zeta_stop = c(-Inf, 0),
            epsilon_go = c(-Inf, 0), epsilon_stop = c(-Inf, 0),
            N = 1,
            df = subset(data, id == i) %>% mutate(id = 1)
          )$data$rt
          
        ) )
      }
      
    ) %>% t()
    
  }
)
