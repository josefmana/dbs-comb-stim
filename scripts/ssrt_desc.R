# This script extracts classical SSRT parameters for description of each participant.
# The script is based on "analyse_stop.R" from Verbruggen.

rm( list = ls() ) # clear the environment

# list required packages into a character object
pkgs <- c( "here", "dplyr", "tidyverse" )

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# prepare folders for results
sapply( c("figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )


# step no.1: prepare the data ----

# read the data
d0 <- read.csv( here("_data","ssrt_raw.csv"), sep = "," )
#d1 <- read.csv( here("_data","ssrt_lab.csv"), sep = "," )
d0 <- d0 %>% mutate( subject = paste(id,sess,sep="_") ) # scaffolding

# exclude practice block
d0 <- subset (d0, block > 0)
#d1 <- subset (d1, block > 0)

# prepare some variables
d0 <- #d1 <-
  
  d0 %>% #d1 %>%
  
  mutate(
    
    # make factors of a couple of variables that may be used later on
    signal = factor( signal, levels = 0:1, labels = c("nosignal","signal"), ordered = T ),
    #signal = factor( signal, levels = c("nosignal","signal"), ordered = T ),
    
    # create new variables for calculation of p(correct)
    acc = ifelse( correct == 4, 1, 0 ), #acc = ifelse( correct == "correct", 1, 0 ),
    miss = ifelse( resp1 == 0, 1, 0 ), #miss = ifelse( is.na(resp1), 1, 0 ),
    presp = ifelse( resp1 > 0, 1, 0 ) #presp = ifelse( !is.na(resp1), 1, 0 ) # means p(response)?
    
  )


# step no.2: do some basic design checks & basic performance ----

# check design
with( d0, table( id, signal, sess ) ) # signal in 1/4 of trials
with( d0, table( id, stim, sess ) ) # the same number of rights and lefts
#with( d1, table( id, signal, sess ) )
#with( d1, table( id, stim, sess ) )

# check overall error rates of all subjects on both no-signal & signal trials
# choosing error rates over accuracy rates for two reasons:
# 1) error rates ought to be smaller and 

# prepare a table of subjects and counterbalancing conditions
t1 <-
  
  # start by preparing an id/counterbalancing part of the table
  with( d1, table( sess, cond, id ) ) %>% # extract frequency table of session/condition combinations for each patient
  as.data.frame() %>%
  filter( Freq > 0 & cond == "exp" ) %>% # keep only rows with experimental condition being presented
  mutate( cbal = ifelse( sess == 1, "+0", "0+" ) ) %>% # prepare counterbalancing labels
  select(id,cbal) %>%
  
  # continue by extracting accuracy rates for control and experimental blocks in nosignal (go) / signal (stop) trials
  left_join(
    
    lapply(
      
      c("nosignal","signal"), # loop through nosignal/signal trials
      function(i)
        
        # extract tables with error rates rounded to nearest percentage point
        # since we want condition-specific row percentages, we multiply by 2 (conditions) * 100 (max percentage)
        with( subset( d1, signal == i ), round( 200 * prop.table( table( id, acc, cond ), margin = 1 ), 0 ) ) %>%
        as.data.frame() %>%
        filter( acc == 0 ) %>% # keep errors (acc == 0) only
        pivot_wider( id_cols = id, names_from = cond, values_from = Freq, names_prefix = paste0(i,"_") )
      
    ) %>%
      
      do.call( left_join, . ) %>%
      relocate( nosignal_exp, .before = signal_exp ),

  )



# check overall accuracy of all subjects on no-signal trials
with( subset( d0, signal == "nosignal" ), round( 200 * prop.table( table( id, acc, sess ), margin = 1 ), 0 ) )
#with( subset( d1, signal == "nosignal" ), round( 200 * prop.table( table( id, acc, sess ), margin = 1 ), 0 ) )
#with( subset( d1, signal == "nosignal" ), round( 200 * prop.table( table( id, acc, cond ), margin = 1 ), 0 ) )

# check overall accuracy of all subjects on signal trials
with( subset( d0, signal == "signal" ), round( 200 * prop.table( table( id, acc, sess ), margin = 1 ), 0 ) )
#with( subset( d1, signal == "signal" ), round( 200 * prop.table( table( id, acc, sess ), margin = 1 ), 0 ) )
#with( subset( d1, signal == "signal" ), round( 200 * prop.table( table( id, acc, cond ), margin = 1 ), 0 ) )
