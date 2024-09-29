# In this script we (i) build a heuristic causal model of system under investigation (via DAG(s)), (ii) build a statistical model based
# on assumed data-generating process and implied dependencies from (i), and (iii) validate the model from (ii) by recovering generating
# parameters when noise is added


# list packages to be used
pkgs <- c( "here", "dplyr", "tidyverse", "dagitty", "ggdag", "ggplot2", "patchwork" )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders for figures and tables to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )


# HEURISTIC CAUSAL MODEL(S) ----

# prepare a plausible DAG (directed acyclic graph)
dag <-
  
  # write the model in dagitty
  dagitty( "dag { motor <- freq -> inhib <- motor <- ldopa -> inhib <- cloc -> motor }" ) %>%
  
  # set-up coordinates  
  `coordinates<-`(
    c(
      list(
        x = c(motor = 0, inhib = 1, freq = 1, cloc = 0, ldopa = 0.5),
        y = c(motor = 1, inhib = 1, freq = 0, cloc = 0, ldopa = 2)
      )
    )
  )

# plot the DAG
f.dag <-
  
  # prepare a list that will contain all the DAGs
  list(
    
    # plot the basic DAG with no backdoors being analysed
    full =
      
      dag %>%
      ggdag() +
      theme_dag() +
      labs( title = "base DAG", subtitle = "(assumed causal structure)" ) +
      theme( plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) )
    
  )

# add the minimal adjustment sets for total and direct effects of stimulation frequency ("freq") combined with
# the contact location ("cloc")
for ( i in c("total","direct") ) {
  
  f.dag[[i]] <-
    
    dag %>%
    ggdag_adjustment_set( exposure = c("cloc","freq"), outcome = "inhib", effect = i, type = "minimal", shadow = F ) +
    theme_dag() +
    theme( legend.position = "none" ) +
    labs( title = paste0(i," effect"), subtitle = "inhib ~ f(cloc, freq)" ) +
    theme( plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) )

}

# arrange the DAGs
with( f.dag, ( full | total | direct ) )

# save it
ggsave( "figs/dag.jpg", dpi = 300, width = 10.5, height = 6.21 )


# DATA-GENERATING MDOEL ----




# SESSION INFO -----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/model_building.txt" )

