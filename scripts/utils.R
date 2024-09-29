# In this script I prepare all in-house functions used thoughout the project.

library(tidyverse)
library(ggh4x)


# ---- SANITY CHECK PLOT ----
sanity_plot <- function(data) data %>%
  
  # some formatting shinanigans
  filter( complete.cases(rt) ) %>%
  mutate( subject = factor(paste0("Subject #",id), levels = unique( paste0("Subject #",id) ), ordered = T) ) %>%
  
  # plot it
  ggplot() +
  aes( x = rt, linetype = as.character(signal) ) +
  geom_density(linewidth = .68) +
  facet_wrap(~ subject, ncol = 3) +
  labs(
    title = "Observed data distributions",
    subtitle = "GO response times (solid) should be generally slower than STOP-RESPOND response times (dashed).",
    x = "Response time (s)",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5)
  )


# ---- PARAMETERS RECOVERY ----
reco_hist <- function(
    
  data, # parameter estimates to be evaluated against a ground truth
  col = c(go = "red", stop = "blue", hist = "grey68"), # colours denoting true data generating values
  quants = c("mean", "sd"), # quantities of interest, will be printed in the listed order
  tit = "Mean/Scale structure" # plot title
  
) {
  
  # pre-process data
  data <- data %>%
    
    filter( parameter %in% quants ) %>%
    mutate( parameter = factor(parameter, levels = quants, ordered = T) )
  
  # plot it  
  data %>%
    
    ggplot() +
    aes(x = value) +
    geom_histogram(fill = col["hist"]) +
    geom_vline(aes(xintercept = true_value), data = data %>% filter(type == "go"), colour = col["go"], linewidth = 1) +
    geom_vline(aes(xintercept = true_value), data = data %>% filter(type == "stop"), colour = col["stop"], linewidth = 1) +
    geom_vline(aes(xintercept = observed), colour = "black", linewidth = 1, linetype = "dashed") +
    facet_grid2(ID ~ type*parameter, scales = "free", independent = T) +
    labs(
      y = NULL,
      title = tit,
      subtitle = "Histograms represent model posteriors, coloured thick lines represent true data-generating values,\nblack dashed lines represent in-sample values from GO trials."
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = .5),
      plot.subtitle = element_text(hjust = .5),
      axis.text.y = element_blank()
    ) %>%
    return()
  
}


# ---- DENSITY OVERLAY OF POSTERIOR PREDICTIONS ----
ppc_density <- function(data, preds, cols, ncols = 4, ndrws = 50) lapply(
  
  1:length(preds),
  function(i)
    
    cbind.data.frame(
      subset(data, id == i),
      preds[[i]][ sample(1:nrow(preds[[i]]), ndrws), subset(data, id == i)$trial ] %>% t()
    )
  
) %>%
  
  # prepare data
  do.call( rbind.data.frame, . ) %>%
  pivot_longer( cols = c( "rt",as.character(1:ndrws) ), names_to = "sample", values_to = "Response time (s)" ) %>%
  mutate(
    source = if_else(sample == "rt", "observed", "predicted"), # source of data (model vs sampling)
    type = paste0(source,"_",signal), # helping variable denoting all types of lines to be printed
    sample = paste0(signal,"_",sample), # grouping variable
    id = factor(paste0("Subject #",id), levels = unique( paste0("Subject #",id) ), ordered = T) # ordered ID
  ) %>%
  
  ggplot() +
  aes(x = `Response time (s)`, colour = type, size = source, group = sample) +
  geom_density() +
  scale_size_manual( values = c(.8,.1) ) +
  scale_colour_manual(values = cols) +
  facet_wrap( ~ id, ncol = ncols, scales = "free" ) +
  labs(
    title = "Posterior predictive check of response time distributions",
    subtitle = "Thick lines represent observed data, thin lines represent model posterior predictions,\nGO trials are presented in red, STOP RESPOND trials in blue"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
