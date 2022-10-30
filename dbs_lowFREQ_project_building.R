# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c( "dplyr", "tidyverse" ) # for data wrangling

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}


# ---- counterbalancing ----

# list all counterbalancing conditions w.r.t. the exposure (i.e., frequency manipulation)
s1 <- c("base/low", "low/base") # the first session: base (high-frequency dorsal stimulation only) vs low-frequency stimulation
s2 <- c("base/low", "low/base") # the second session: the same conditions as the first session
int <- c("low", "base") # the interim, either combined low-ventral and high-dorsal or high-frequency (base) stimulation only

# generate all counterbalancing conditions across sessions and measurements
cb.long <- expand.grid(s1,int,s2) %>%
  `colnames<-`( c("session_1", "interim", "session_2") ) %>%
  mutate( cb_long = paste0( "cond_", 1:nrow(.) ),
          change_freq = ifelse( session_1 == session_2, 0, 1)
          )

# list all counterbalancing conditions for the cognitive assessments
block <- c("srt/vf", "vf/srt")

