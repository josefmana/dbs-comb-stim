# In this script we (i) build a heuristic causal model of system under investigation (via DAG(s)), (ii) build a statistical model based
# on assumed data-generating process and implied dependencies from (i), and (iii) validate the model from (ii) by recovering generating
# parameters when noise is added


# list packages to be used
pkgs <- c( "rstudioapi", "dplyr", "tidyverse", "dagitty", "ggdag", "ggplot2", "patchwork" )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# create folders "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("figs", "tabs", "sess"), function(i) if( !dir.exists(i) ) dir.create(i) )


# ---- HEURISTIC CAUSAL MODELS (DAGs) ----

# prepare a plausible DAG
dag <- dagitty( "dag { motor <- freq -> inhib <- motor <- ldopa -> inhib <- cloc -> motor }" ) %>%
  
  # set-up coordinates  
  `coordinates<-`( c( list( x = c(motor = 0, inhib = 1, freq = 1, cloc = 0, ldopa = 0.5),
                            y = c(motor = 1, inhib = 1, freq = 0, cloc = 0, ldopa = 2)
                            ) ) )

# plot the DAG
f.dag <- list( full = dag %>% ggdag() + theme_dag() +
                 labs( title = "base DAG", subtitle = "(assumed causal structure)" ) +
                 theme( plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) ) )

# add the minimal adjustment sets for total and direct effects of stimulation frequency ("freq") combined with
# the contact location ("cloc")
for ( i in c("total","direct") ) f.dag[[i]] <- dag %>%
  ggdag_adjustment_set( exposure = c("cloc","freq"), outcome = "inhib", effect = i, type = "minimal", shadow = F ) +
  theme_dag() + theme( legend.position = "none" ) +
  labs( title = paste0(i," effect"), subtitle = "inhib ~ f(cloc, freq)" ) +
  theme( plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) )

# arrange the DAGs
( f.dag$full | f.dag$total | f.dag$direct )

# save it
ggsave( "figs/dag.jpg", dpi = 300, width = 1.1*10.5, height = 1.1*6.21 )


# ---- DATA-GENERATING MDOEL ----




# ---- SESSION INFO -----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/model_building.txt" )

