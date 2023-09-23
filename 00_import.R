# This script extracts and formats data from raw tables to analysis-ready outcomes.

# list required packages into a character object
pkgs <- c( "here", "dplyr", "tidyverse" )

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# prepare a folder for sessions info and pre-processed data
sapply( c("sess","data") , function(i) if( !dir.exists(i) ) dir.create(i) )

# extract important headers
paths <- here( "raw", list.files( here("raw") ) ) # paths to each single patient/session data set
sid <- paths %>% gsub("\\..*","",.) %>% gsub(".*-","",.) %>% sub("_.*", "", . ) %>% unique() # subject IDs

# read and wrangle the data
d0 <- lapply(
  
  # loop through all subjects
  setNames(sid,sid), function(i)
    
    lapply(
      # loop through sessions
      1:2, function(j)
        # read and wrangle the data
        read.table( paths[ grepl( i, paths ) & grepl( paste0(j,"."), paths, fixed=T ) ], header = T ) %>%
        mutate( sess = j, .before = 1 ) %>% # add session
        mutate( id = i, .before = 1 ) # add subject ID
    ) %>%
    # add data from both sessions to a single data frame
    do.call( rbind.data.frame, . )

) %>%
  # put together data sets of all subjects
  do.call( rbind.data.frame, . ) %>%
  `rownames<-`( 1:nrow(.) )
