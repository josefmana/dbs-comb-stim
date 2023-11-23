# This script extracts and formats data from raw tables to analysis-ready outcomes.

# list required packages into a character object
pkgs <- c( "here", "readxl", "janitor", "dplyr", "tidyverse" )

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# prepare a folder for sessions info and pre-processed data
if( !dir.exists("_data") ) dir.create("_data")


# DATA PATHWAYS ----

# paths to each single patient/session data set
paths <- here( "_raw", "ssrt", list.files( here("_raw","ssrt") ) )

# mapping of subject/session pairs to their respective paths
d0 <-
  paths %>%
  sub( "\\..*", "", . ) %>% # drop file type ending
  sub( ".*-", "", . ) %>% # drop the path info up to "stop-"
  tolower() %>% # put all characters to lowercase
  strsplit("_") %>% # split variable type (id, surname, forename, session)
  do.call( rbind.data.frame, . ) %>% # create a single data set
  `colnames<-`( c("id","surname","forename","session") ) %>% # label the columns appropriately
  cbind.data.frame( . , path = paths ) %>% # add path information to each patient/session combination
  mutate( session = as.integer( gsub("\\D","",session) ) ) # change session number to a number


# COUNTERBALANCING ----

# read the table
cb <-
  read_xlsx( here( "_raw","rand","counterkey.xlsx") ) %>% # read the data
  
  # renaming variables
  rename( "session1" = "1.vyšetření" ) %>%
  rename( "session2" = "2.vyšetření" ) %>%
  rename( "surname" = "jméno" ) %>%
  
  # transforming variables
  mutate( across( contains("session"), ~ ifelse( is.na(.x), 0, 1 ) ) ) %>%
  mutate( surname = sub( "_.*", "", make_clean_names(surname) ) )

# check all surnames are included in both tables
all( cb$surname %in% d0$surname ) # TRUE

# NOTE: THIS PORTION OF THE CODE WILL NEED SOME BETTER CHECKING BECAUSE IT WOULD NOT WORK IF THERE WERE TWO PATIENTS
# WITH THE SAME SURNAME IF USED IN ITS CURRENT FORM

# add indicator of experimental session to the pathway table
d0$cfs <- sapply( 1:nrow(d0), function(i) with( cb, get( paste0("session",d0$session[i] ) )[surname == d0$surname[i] ] ) )


# SSRT DATA ----

# ---- raw in-numbers data ----

# read and wrangle the data
d1 <-
  
  with(
    
    d0,
    lapply(
      
      # loop through all subjects
      setNames( unique(id), unique(id) ),
      function(i)
        
        lapply(
          
          # loop through sessions
          1:2,
          function(j)
          
            read.table( path[ id == i & session == j ], header = T ) %>% # read the data
            mutate( exp = cfs[ id == i & session == j ] , .before = 1 ) %>% # add experimental condition
            mutate( sess = j, .before = 1 ) %>% # add session
            mutate( id = i, .before = 1 ) # add subject ID
          
        ) %>%
        
        # add data from both sessions to a single data frame
        do.call( rbind.data.frame, . )

    ) %>%
      
      # put together data sets of all subjects
      do.call( rbind.data.frame, . ) %>%
      `rownames<-`( 1:nrow(.) )
      
  )

# save this as a raw numbers-based data set as csv
write.table( x = d1, file = "_data/ssrt_raw.csv", sep = ",", row.names = F, quote = F )


# ---- tidy labelled data ----

# re-code selected columns
d2 <- 
  
  d1 %>%
  
  # add condition label
  mutate( cond = ifelse( exp == 1, "exp", "ctrl" ), .after = exp ) %>%
  
  # transform number variables to labels
  mutate(
    
    stim = case_when( stim == 1 ~ "left", stim == 2 ~ "right" ), # stimulus identity
    signal = ifelse( signal == 1, "signal", "nosignal" ), # stop-signal presence
    reqSOA = ifelse( reqSOA == 0, NA, reqSOA ), # required stimulus onset asynchrony
    
    # correctness of the response
    correct = case_when(
      correct == 4 ~ "correct",
      correct == 3 ~ "signal-respond",
      correct == 2 ~ "incorrect",
      correct == 1 ~ "missed"
    ),
    
    # response recorded
    across(
      starts_with("resp"),
      ~ case_when(
        .x == 0 ~ NA, # miss
        .x == 1 ~ "left",
        .x == 2 ~ "right",
        .x == 3 ~ "down"
      )
    ),
    
    across( starts_with("rt"), ~ ifelse( .x == 0, NA, .x )  ) # response times

  )

# save this as a labelled data set as csv
write.table( x = d2, file = "_data/ssrt_lab.csv", sep = ",", row.names = F, quote = F )
