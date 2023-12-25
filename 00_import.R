# This script extracts and formats data from raw tables to analysis-ready outcomes.

# clear the environment
rm( list = ls() )
gc()

# list required packages into a character object
pkgs <- c( "here", "readxl", "janitor", "dplyr", "tidyverse", "lubridate" )

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
  mutate( surname = sub( "_.*", "", make_clean_names(surname) ) ) %>%
  mutate( datum = ifelse( datum == "25.5.52023", "25.5.2023", datum ) ) # manually change typos in dates

# check all surnames are included in both tables
all( cb$surname %in% d0$surname ) # TRUE

# NOTE: THIS PORTION OF THE CODE WILL NEED SOME BETTER CHECKING BECAUSE IT WOULD NOT WORK IF THERE WERE TWO PATIENTS
# WITH THE SAME SURNAME IF USED IN ITS CURRENT FORM

# add indicator of experimental session and date of the session to the pathway table
d0$cfs <- sapply( 1:nrow(d0), function(i) with( cb, get( paste0("session",d0$session[i] ) )[surname == d0$surname[i] ] ) )
d0$date <- as.Date( sapply( 1:nrow(d0), function(i) with( cb, datum[surname == d0$surname[i] ] ) ), "%d.%m.%Y" )

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


# DEMOGRAPHIC DATA ----

# read and prepare the data
d3 <-
  
  # read it
  read.csv( here("_raw","demo","redcap_data.csv"), sep = "," ) %>%
  
  # keep only data of included patients
  mutate( id = tolower(study_id) ) %>%
  filter( id %in% unique(d0$id) ) %>%
  
  # select variables of interest
  select(
    
    id, dob, sex, # demographic variables
    type_pd, hy_stage, rok_vzniku_pn, asym_park, # PD-specific variables
    surgery_date, datum_stim, # surgery data
    
    # stimulation data
    datum,
    dbs_dex, current_right, duration_right, frequency_right, impedance_right, # right side
    dbs_sin, current_left, duration_left, frequency_left, impedance_left, # left side
    
    # cognition
    
    datum_drs, drsii_total, moca_e04359, nart_7fd846, # level-I neuropsychology
    
    # level-II neuropsychology
    datum_neuropsy_23afdc,
    lns, ds_b, corsi_b, tmt_a, pst_d, # attention & working memory
    tol_anderson, tmt_b, pst_w, pst_c, vf_skp, cf, # executive function
    sim, bnt_60, # language
    ravlt_irs, ravlt_b, ravlt_6, ravlt_30, ravlt_drec50, ravlt_drec15, bvmt_irs, bvmt_30, bvmt_drec, # memory
    jol, clox_i, # visuospatial function
    
    # psychomotor speed/hand-eye coordination
    gp_r, gp_l
    
  ) %>%
  
  # re-code some
  mutate(
    
    # demographics and Parkinson's related variables
    sex = case_when( sex == 0 ~ "female", sex == 1 ~ "male" ),
    type_pd = case_when( type_pd == 1 ~ "tremordominant", type_pd == 2 ~ "akinetic-rigid" ),
    asym_park = case_when( asym_park == 1 ~ "right", asym_park == 2 ~ "left" ),
    
    # dates
    across( contains("dob") | contains("dat"), ~ as.Date(.x) ),
    
    # stimulation parameters 
    across( contains("impedance"), ~ as.numeric( sub( "ok", NA, .x ) ) ),
    across( contains("current"), ~ as.numeric( sub( ",", ".", .x ) ) ),
    across( contains("duration"), ~ as.numeric(.x) ),
    across( contains("frequency"), ~ as.numeric(.x) )

  )
  
# add time differences
for ( i in names(d3)[ grepl( "dob|dat", names(d3) ) ] ) {
  for ( j in 1:nrow(d3) ) {
    
    d3[ j , paste0(i,"_diff") ] <-
      
      time_length( difftime( d0[ with( d0, session == 1 & id == d3$id[j] ), "date" ], d3[ j , i ] ), "years" )
    
  }
}


# SESSION INFO -----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "import_envir.txt" )
