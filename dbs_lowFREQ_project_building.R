# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c( "dplyr", "tidyverse" ) # for data wrangling

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )


# ---- counterbalancing ----

# list all counterbalancing conditions w.r.t. the exposure (i.e., frequency manipulation)
s1 <- c("base/low", "low/base") # the first session: base (high-frequency dorsal stimulation only) vs low-frequency stimulation
s2 <- c("base/low", "low/base") # the second session: the same conditions as the first session
int <- c("low", "base") # the interim, either combined low-ventral and high-dorsal or high-frequency (base) stimulation only

# generate all counterbalancing conditions across sessions and measurements
cb.long <- expand.grid(s1,int,s2) %>%
  `colnames<-`( c("session_1", "interim", "session_2") ) %>%
  write.table( . , file = "tables/counterbalancing.csv", sep = ",", row.names = F )


# ----------- session info -----------

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/dbs_lowFREQ_project_building.txt" )

