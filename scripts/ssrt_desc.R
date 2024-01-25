# This script extracts classical SSRT parameters for description of each participant.
# The script is based on "analyse_stop.R" from Verbruggen.

rm( list = ls() ) # clear the environment

# list required packages into a character object
pkgs <- c("here","tidyverse","ggplot2","patchwork")

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# prepare folders for results
sapply( c("figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )


# step no.1: prepare the data ----

# read the data
d1 <- read.csv( here("_data","ssrt_lab.csv"), sep = "," )

# exclude practice block
d1 <- subset (d1, block > 0)

# prepare some variables
d1 <- #d0 <-
  
  d1 %>% #d0 %>%
  
  mutate(
    
    # make factors of a couple of variables that may be used later on
    signal = factor( signal, levels = c("nosignal","signal"), ordered = T ),
    
    # create new variables for calculation of p(correct)
    acc = ifelse( correct == "correct", 1, 0 ),
    miss = ifelse( is.na(resp1), 1, 0 ),
    presp = ifelse( !is.na(resp1), 1, 0 ) # means p(response)?
    
  )

# prepare a table with id/counterbalancing pairs
t0 <-
  
  with( d1, table( sess, cond, id ) ) %>% # extract frequency table of session/condition combinations for each patient
  as.data.frame() %>%
  filter( Freq > 0 & cond == "exp" ) %>% # keep only rows with experimental condition being presented
  mutate( cbal = ifelse( sess == 1, "+0", "0+" ) ) %>% # prepare counterbalancing labels
  select(id,cbal)

# step no.2: do some basic design checks & basic performance ----

# check design
with( d1, table( id, signal, sess ) ) # signal in 1/4 of trials
with( d1, table( id, stim, sess ) ) # the same number of rights and lefts

# check overall error rates of all subjects on both no-signal & signal trials
# choosing error rates over accuracy rates for two reasons:
# 1) error rates ought to be smaller (easier to spot high values) and 
# 2) error rates has subtypes that can be analysed in next sections

# prepare a table of subjects and counterbalancing conditions
t1 <-
  
  # extract accuracy rates for control and experimental blocks in nosignal (go) / signal (stop) trials
  t0 %>%
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

# save it
write.table( t1, here("tabs","error_tab.csv"), sep = ",", row.names = F, quote = F )


# step no. 3: analyse no-signal data ----

# prepare a table summarising accuracy, miss rate and response times of no-signal (go) trials
t2 <-
  
  t0 %>%
  left_join(
    
    d1 %>%
      
      filter( signal == "nosignal" ) %>% # only nosignal/go trials considered
      group_by( id, cond ) %>% # summarise for each patient/condition combination
      summarise(

        `p(correct)` = mean( acc[presp == 1] ), # accuracy, i.e. p(correct) (omitting trials without a response or anticpatory responses )
        `p(miss)` = mean( miss ), # miss rate
        
        # response times for correct responses
        RT =
          paste0(
            round( mean( rt1[correct == "correct"], na.rm = T ), 0 ), " (",
            round( sd( rt1[correct == "correct"], na.rm = T ), 0 ), ")"
          )

      ) %>%
      
      # tidy it up
      pivot_wider( values_from = c(`p(correct)`, `p(miss)`, RT ), names_from = cond ) %>%
      relocate( "p(miss)_ctrl", .after = "p(correct)_ctrl" ) %>%
      relocate( "RT_ctrl", .after = "p(miss)_ctrl" )
    
  )

# save it
write.table( t2, here("tabs","nosignal_tab.csv"), sep = ",", row.names = F, quote = F )

# prepare a graphical representation of t2
f1 <-
  
  lapply(
    
    setNames( c("correct","miss"), c("correct","miss") ),
    function(i)
      
      t2 %>%
      
      # prepare the data
      select( id, contains(i) ) %>%
      pivot_longer( -id, names_to = "Condition:", values_to = paste0("Pr(",i,")") ) %>%
      rename( "ID" = "id" ) %>%
      mutate(
        `Condition:` =
          case_when( grepl("exp",`Condition:`) ~ "experimental", grepl("ctrl",`Condition:`) ~ "control" ) %>%
          factor( levels = c("control","experimental"), ordered = T )
      ) %>%
      
      # plotting proper
      ggplot() +
      aes( x = ID, y = get( paste0("Pr(",i,")") ), fill = `Condition:` ) +
      geom_bar( stat = "identity", position = position_dodge(), width = .5 ) +
      labs( x = NULL, y = paste0("Pr(",i,")") ) +
      coord_cartesian( ylim = case_when( i == "correct" ~ c(.9,1), .default = NULL ) ) +
      theme_bw( base_size = 14 ) +
      theme( legend.position = "none" ) + # keep legend of the response time plot only
      scale_fill_manual( values = cbPal[c(1,2)] )

  )

# add error bar plots for response times
f1$resp <-
  
  t2 %>%
  
  # prepare the data
  select( id, contains("RT") ) %>%
  pivot_longer( -id, names_to = "Condition:", values_to = "Response time (ms)" ) %>%
  rename( "ID" = "id" ) %>%
  mutate(
    `Condition:` =
      case_when( grepl("exp",`Condition:`) ~ "experimental", grepl("ctrl",`Condition:`) ~ "control" ) %>%
      factor( levels = c("control","experimental"), ordered = T ),
    sd = as.numeric( gsub( "\\D", "", sub( ".* ", "", `Response time (ms)` ) ) ),
    `Response time (ms)` = as.numeric( sub( " .*", "", `Response time (ms)`) ),
    low = `Response time (ms)` - sd,
    upp = `Response time (ms)` + sd
  ) %>%
  
  # plotting proper
  ggplot() +
  aes( x = ID, y = `Response time (ms)`, colour = `Condition:` ) +
  geom_point( size = 4, position = position_dodge( width = 1/3) ) +
  geom_linerange( aes( ymin = low, ymax = upp ), position = position_dodge( width = 1/3), size = 4, alpha = .5 ) +
  theme_bw( base_size = 14 ) + theme( legend.position = "bottom" ) +
  scale_colour_manual( values = cbPal[c(1,2)] )

# put them to a single figure
with( f1, correct/miss/resp ) +  plot_annotation( tag_levels = "a" )

# save the figure
ggsave( here("figs","nosignal_fig.jpeg"), width = 8.5, height = 11 )


# step no. 4: analyse signal data ----

# prepare a function to calculate all signal data at once
funcSignal <-
  
  function(d) {
    
    # calculate prespond & mean stop-signal delay
    with(
      subset( d, signal == "signal" ), { # use signal only data
        presp_m <<- mean(presp)
        ssd_m <<- mean(trueSOA)
        ssd_v <<- sd(trueSOA)
      }
    )
    
    # nth response time & mean (sd) go response time
    with(
      subset( d, signal == "nosignal" & !is.na(resp1) ), { # extract no-signal data (drop trials with missed responses)
        nthRT <<- quantile( rt1, probs = presp_m, type = 6, na.rm = T, names = F )
        goRT_m <<- mean( rt1, na.rm = T )
        goRT_v <<- sd( rt1, na.rm = T )
      }
    )
    
    # SSRT = nthRT - ssd
    ssrt <- nthRT - ssd_m
    
    # signal-respond RT
    with(
      subset( d, signal == "signal" & presp == 1), {
        sRT_m <<- mean( rt1, na.rm = T )
        sRT_v <<- sd( rt1, na.rm = T )
      }
    )
    
    # return
    return(
      c(
        ssd = paste0( round(ssd_m,0), " (", round(ssd_v,0), ")" ), # Stop-Signal Delay
        nthRT = round( nthRT, 0 ), # nth go reaction time
        ssrt = round( ssrt, 0 ), # Stop-Signal Response Time
        sRT = paste0( round(sRT_m,0), " (", round(sRT_v,0), ")" ), # signal response time
        raceTest = round( goRT_m - sRT_m, 0 ) # difference between mean go vs stop response times
      )
      
    )  

  }

# pull it all into a single nice neat table
t3 <-
  
  t0 %>%
  left_join(
    
    lapply(
      
      unique(d1$cond), # loop through experimental conditions
      function(i)
        
        sapply(
          
          unique(d1$id), # loop through patients
          function(j)
            
            # bind patient's ID to their signal data summaries
            c( id = j, funcSignal( d1[ with( d1, id == j & cond == i ) , ] ) )
          
        ) %>%
        
        # tidy it up
        t() %>%
        as.data.frame() %>%
        mutate( cond = i, .after = 1 )
      
    ) %>%
      
      # pull control and experimental data to a single file
      do.call( rbind.data.frame, . ) %>%
      pivot_wider( values_from = -c(id,cond), names_from = cond, id_cols = id ),
    
    by = "id"
    
  ) %>%
  
  # change column order such that experimental conditions are grouped together
  relocate( "ssd_exp", .after = "ssrt_ctrl" ) %>%
  relocate( "nthRT_exp", .after = "ssd_exp") %>%
  relocate( "raceTest_ctrl", .before = "ssd_exp" )

# save it
write.table( t3, here("tabs","signal_tab.csv"), sep = ",", row.names = F, quote = F )

# NEED TO DO: VISUALISE SIGNAL TRIALS



# SESSION INFO -----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = here("scripts","ssrt_desc_envir.txt") )
