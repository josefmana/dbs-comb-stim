# This file is for easier data writing of metadata for the dbs_combSTIM project

# list packages to be used
pkgs <- c("rstudioapi", "dplyr", "tidyverse")

# load or install each of the packages as needed
for (i in pkgs) {
  if (i %in% rownames (installed.packages()) == F) install.packages(i) # install if it ain't installed yet
  if ( i %in% names (sessionInfo() $otherPkgs) == F ) library(i, character.only = T ) # load if it ain't loaded yet
}

# set working directory (works in RStudio only)
setwd( paste0( dirname(getSourceEditorContext()$path), "/_no_github/data" ) )

# read current data sets
speech <- read.csv( "dbs_combSTIM_speech_metadata.csv", sep = ";" )
stai <- read.csv( "dbs_combSTIM_stai.csv", sep = "," )


# ---- list patients id and name ----
id = NA
name = NA


# ---- add the patient to STAI ----

# prepare the data to be added
stai_start = rep(NA,2) # the time STAI started (a vector of two)
stai_score = NA # patient STAI score (collapsed over all sessions)

# add new rows to the stai data frame
stai <- stai %>%
  # add the new patient
  rbind.data.frame( 
    # prepare a data frame for the new patient
    data.frame( id = id, name = name, # identificators
                session = c( rep(1,40) , rep(2,40) ), start = c( rep( stai_start[1], 40), rep( stai_start[2], 40) ), # session info
                scale = rep( c( rep("X1",20), rep("X2",20 ) ), 2 ), item = rep( paste0("i",1:20), 4 ), # scale info
                score = stai_score # scores
    )
  )


# --- add the patient's speech metadata ----

# prepare the data to be added
speech_cb = rep(NA,4) # the order in which verbal fluency stimuli were presented (a vector of four)
speech_start = rep(NA,8) # starting time of each task in the speech blocks

# add new rows to the speech data frame
speech <- speech %>%
  # add a data.frame for new patient below the existing ones
  rbind.data.frame(
    # prepare a data frame for the new patient
    data.frame( id = id, name = name, # identifiers
                session = c( rep(1,4), rep(2,4) ), # session info
                task = rep( c("letter_fluency","category_fluency","reading","speaking"), 2 ), task_order = rep( 1:4, 2 ),
                task_spec = c( speech_cb[1:2], "seed", "cinderella", speech_cb[3:4], "travel", "gingerbread_house"), # task specification
                file = rep( paste0( id, "_", name, "_session", 1:2, ".wav" ), 4) %>% sort(),
                start = speech_start
    )
  )


# ---- check & save the new files ----

for( i in c("stai","speech") ) View( get(i) )

# if A-OK, continue by saving as .csv
for ( i in c("stai","speech") ) write.table( get(i),
                                             paste0( "dbs_combSTIM_", ifelse( i == "stai", i, paste0(i, "_metadata") ), ".csv" ),
                                             sep = ifelse( i == "stai", ",", ";" ), row.names = F, quote = F
                                             )
