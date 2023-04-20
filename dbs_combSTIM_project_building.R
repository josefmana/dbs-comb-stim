# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c( "dplyr", "tidyverse", # for data wrangling
           "dagitty", "ggdag", # heuristic causal models
           "ggplot2", "patchwork" # plotting
           )

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("figs", "tabs", "sess"), function(i) if( !dir.exists(i) ) dir.create(i) )


# ---- counterbalancing ----

# list all counterbalancing conditions w.r.t. the exposure (i.e., frequency manipulation)
s1 <- c("dorsal/combined", "combined/dorsal") # the first session: base (high-frequency dorsal stimulation only) vs combined stimulation
s2 <- c("dorsal/combined", "combined/dorsal") # the second session: the same conditions as the first session
int <- c("combined", "dorsal") # the interim, either combined low-ventral and high-dorsal or high-frequency (base) stimulation only

# generate all counterbalancing conditions across sessions and measurements
cb.long <- expand.grid( s1, int, s2 ) %>% `colnames<-`( c("session_1", "interim", "session_2") )

# save as .csv
write.table( cb.long , file = "tabs/counterbalancing.csv", sep = ",", row.names = F )


# ---- trajectories ----

# based on Filip Ruzicka's e-mail from 2023-02-11 (18:39), working on a more complex counterbalancing scheme
# start by listing all possible combinations for the first session across left/right and base/comb stimulation
# this one assumes that the patient will be in the interim stimulate by the last mode used in session_1
traj <- expand.grid( side = c("left","right"),
                     session_1 = c("dors-comb","comb-dors"),
                     session_2 = c("dors-comb","comb-dors")
                     ) %>%
  mutate( interim = sub( ".*-", "", session_1 ), .before = "session_2" )

# save as .csv
write.table( traj , file = "tabs/trajectories.csv", sep = ",", row.names = F )

# try another way
expand.grid( `1a` = c("left","right","base"),
             `1b/2a` = c("left","right","base"),
             `2b/3a` = c("left","right","base"),
             `3b/4a` = c("left","right","base")
             )  %>%
  mutate( Var5 = ifelse( `1b/2a` == `2b/3a` | `1b/2a` == `3b/4a` | `2b/3a` == `3b/4a` , NA, 1 ) ) %>%
  na.omit() %>% `rownames<-`( 1:nrow(.) ) %>% select(-5)

# ---- heuristic causal models (DAGs) ----

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


# ----------- session info -----------

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/dbs_combSTIM_project_building.txt" )

