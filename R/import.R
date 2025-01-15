#
# This script extracts and formats patient data from raw tables to analysis-ready outcomes.
#

#
# LIST PATHWAYS TO DATA ----
data_paths <- function(folder) {
  
  paths <- here( "_raw", folder, list.files( here("_raw", folder) ) ) # extract raw paths
  
  # pre-format
  paths %>%
    sub("\\..*", "", . ) %>%                                     # drop file type ending
    sub(".*-"  , "", . ) %>%                                     # drop the path info up to "stop-"
    tolower() %>%                                                # put all characters to lowercase
    strsplit("_") %>%                                            # split variable type (id, surname, forename, session)
    do.call( rbind.data.frame, . ) %>%                           # create a single data set
    `colnames<-`( c("id","surname","forename","session") ) %>%   # label the columns appropriately
    cbind.data.frame( . , path = paths ) %>%                     # add path information to each patient/session combination
    mutate( session = as.integer( gsub("\\D","",session) ) ) %>% # change session number to a number
    mutate( id = sub( "ipn", "IPN", id ) ) %>%                   # return ids to uppercase
    return()
  
}

#
# EXTRACT COUNTERBALANCING ----
extract_counterbalancing <- function(.file) .file %>%
  
  mutate( # re-code sessions (+LFS -> 1, otherwise -> 0)
    across(
      contains("round"),
      ~ ifelse( is.na(.x), 0, 1 )
    )
  ) %>%
  mutate( # get rid of diacritics
    surname = sub(
      pattern = "_.*",
      replacement = "",
      x = make_clean_names(surname)
    )
  ) %>%   
  mutate( # manually change typos in dates
    date = if_else(
      condition = date == "25.5.52023",
      true = "25.5.2023",
      false = date
    )
  )

#
# IMPOORT SSRT DATA ----
get_data <- function(.paths, .cb, .demofile, type = "raw") {
  
  # check all surnames are included in both tables
  if(all( .cb$surname %in% .paths$surname ) == F) return("Surnames do not match between data pahts and counterbalancing file!")
  
  # if everything is all-right, read the data
  else {
    
    # extract patient ids
    pats <- unique(.paths$id)
    
    # add indicator of experimental session and date of the session to the pathway table
    .paths <- .paths %>%
      
      mutate(
        cfs = unlist(
          sapply(
            X = 1:nrow(.),
            FUN = function(i) .cb[.cb$surname == surname[i], paste0("round_", session[i]) ]
          ),
          use.names = F
        ),
        date = as.Date(
          x = unlist(
            x = sapply(
              X = 1:nrow(.),
              FUN = function(i) .cb[.cb$surname == surname[i], "date"]
            ),
            use.names = F
          ),
          format = "%d.%m.%Y"
        )
      )
    
    ## ---- raw in-numbers data ----
    raw <- with(
      
      .paths, lapply(
        
        # loop through all subjects
        X = set_names( unique(id) ),
        FUN = function(i) lapply(
          
          # loop through sessions
          X = 1:2,
          FUN = function(j)

            read.table(path[ id == i & session == j ], header = T) %>%    # read the data
            mutate(exp = cfs[ id == i & session == j ] , .before = 1) %>% # add experimental condition
            mutate(sess = j, .before = 1) %>%                             # add session
            mutate(id = i, .before = 1)                                   # add subject ID

        ) %>%

          # add data from both sessions to a single data frame
          do.call( rbind.data.frame, . )

      )  %>%

        # put together data sets of all subjects
        do.call( rbind.data.frame, . ) %>%
        `rownames<-`( 1:nrow(.) )

    )
  
    ## ---- tidy labelled data ----
    labelled <- raw %>%
      
      # add condition label
      mutate(cond = if_else(exp == 1, "exp", "ctrl"), .after = exp) %>%
      
      # transform number variables to labels
      mutate(
        
        stim = case_when(stim == 1 ~ "left", stim == 2 ~ "right"), # stimulus identity
        signal = ifelse(signal == 1, "signal", "nosignal"),        # stop-signal presence
        reqSOA = ifelse(reqSOA == 0, NA, reqSOA),                  # required stimulus onset asynchrony

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
        
        # response times
        across( starts_with("rt"), ~ if_else(.x == 0, NA, .x)  )
        
      )
    
    ## ---- demographic data ----
    demography <- read.csv(.demofile, sep = "," ) %>%
      
      # keep only data of included patients
      filter(study_id %in% pats) %>%
      
      # select variables of interest
      select(
        
        study_id, redcap_event_name, dob, sex,       # demographic variables
        type_pd, hy_stage, rok_vzniku_pn, asym_park, # PD-specific variables
        surgery_date, datum_stim,                    # surgery data
        
        # stimulation data
        datum,
        current_right, duration_right, frequency_right, impedance_right, # right side
        current_left, duration_left, frequency_left, impedance_left,     # left side
        
        # cognition
        datum_drs, drsii_total, moca_e04359, nart_7fd846, # level-I neuropsychology
        
        # level-II neuropsychology
        datum_neuropsy_23afdc,
        lns, ds_b, corsi_b, tmt_a, pst_d,                                                                # attention & working memory
        tol_anderson, tmt_b, pst_w, pst_c, vf_skp, cf,                                                   # executive function
        sim, bnt_60,                                                                                     # language
        ravlt_irs, ravlt_b, ravlt_6, ravlt_30, ravlt_drec50, ravlt_drec15, bvmt_irs, bvmt_30, bvmt_drec, # memory
        jol, clox_i,                                                                                     # visuospatial function
        
        # psychomotor speed/hand-eye coordination
        gp_r, gp_l
        
      ) %>%
      
      # rename some variables
      rename(
        
        # demographics/helpers
        "id" = "study_id",
        "event" = "redcap_event_name",
        "pd_dur" = "rok_vzniku_pn",
        
        # date variables
        "age_years" = "dob",
        "elsurg_years" = "surgery_date",
        "stimsurg_years" = "datum_stim",
        "stim_years" = "datum",
        "drs_years" = "datum_drs",
        "neuropsy_years" = "datum_neuropsy_23afdc",
        
        # neuropsychology
        "moca" = "moca_e04359",
        "nart" = "nart_7fd846",
        "tol" = "tol_anderson"
        
      ) %>%
      
      # re-code some
      mutate(
        
        # demographics and Parkinson's related variables
        event     = sub("_arm_1", "", event) %>% sub("nvtva_", "", .) %>% sub("operace", "surgery", .),
        pd_dur    = 2023 - pd_dur,
        sex       = case_when(sex == 0 ~ "female", sex == 1 ~ "male"),
        type_pd   = case_when(type_pd == 1 ~ "tremor-dominant", type_pd == 2 ~ "akinetic-rigid"),
        asym_park = case_when(asym_park == 1 ~ "right", asym_park == 2 ~ "left"),
        
        # stimulation parameters 
        across( contains("impedance"), ~ as.numeric( sub("ok", NA, .x) ) ),
        across( contains("current")  , ~ as.numeric( sub(",", ".", .x) ) ),
        across( contains("duration") , ~ as.numeric(.x) ),
        across( contains("frequency"), ~ as.numeric(.x) )
        
      ) %>%
      
      # drop refresh assessments
      filter( event != "refresh" )
    
    # add time differences
    for ( i in which( grepl( "years", names(demography) ) ) ) demography[ , i] <- unlist(
      
      sapply(
        1:nrow(demography),
        function(j) time_length(
          x = difftime(
            time1 = .paths[with( .paths, session == 1 & id == demography$id[j] ), "date"],
            time2 = as.Date( demography[j, i] )
          ),
          unit = "years"
        )
      ),
      use.names = F
      
    )
    
    # move post-surgery variables up
    for(i in pats) {
      
      # move the surgery dates
      demography[ with(demography, id == i & event == "screening"), c("elsurg_years","stimsurg_years") ] <- 
        demography[ with(demography, id == i & event == "surgery"), c("elsurg_years","stimsurg_years") ]
      
      # move the closest stimulation parameters and post-surgery DRS-2
      j <- which( with( demography, id == i & abs(stim_years) == min( abs( stim_years[ id == i ] ), na.rm = T ) ) ) # extract row number of the closest stimulation parameters
      demography[ with( demography, id == i & event == "screening"), 11:19 ] <- demography[ j , 11:19 ]
      demography[ with( demography, id == i & event == "screening"), c("drs_years_post","drs_post") ] <- demography[ j , c("drs_years","drsii_total") ]
      
      # remove post-surgery rows
      demography <- demography[ -which( with( demography, id == i & grepl( "surgery|r1|r3", event ) ) ),  ]
      
    }
    
    # finishing touches
    demography <- demography %>%
      
      select(-event) %>%                          # drop column with event-name
      relocate(drs_years_post, .after = moca) %>% # relocate post-surgery DRS
      relocate(drs_post, .after = drs_years_post)
    
    # return selected file
    return( get(type) )
    
  }
}
