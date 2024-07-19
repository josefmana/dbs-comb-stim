# This script is used to run ExGaussian model of SSRT data implemented in the 'dmc' package.

rm( list = ls() )

# For this script to run, the 'dmc' package from osf.io/pbwx8 need to be put into a '_dcm' folder

library(here)
library(tidyverse)


# DATA PRE-PROCESSING ----

rtmax <- 1.5 # maximum time allow was 1.5 s
d0 <- read.csv( here("_data","ssrt_lab.csv"), sep = "," ) # read data

# pivot it wider
d1 <-
  d0 %>%
  mutate( cens = ifelse( correct == "missed", "right", "none" ), max = rtmax ) %>% # add censoring info and maximum possible RTs
  filter( block != 0 ) %>%
  select( id, block, cond, signal, rt1, trueSOA, cens, max ) %>%
  rename( "rt" = "rt1", "ssd" = "trueSOA" ) %>%
  mutate( across( c("rt","ssd"), ~ .x/1e3 ) ) # re-scale from ms to s

# separate files for "go" and "stop" trials
#Dgo <- d1 %>% filter( signal == "nosignal" )
#Dsr <- d1 %>% filter( signal == "signal" & !is.na(rt) )
#Dna <- d1 %>% filter( signal == "signal" & is.na(rt) ) %>% mutate( rt = rtmax, cens = "right" )

# prepare data sets for model
#GOcon <- subset( Dgo, complete.cases(rt) & cond == "ctrl" ) %>% mutate( id = as.integer( as.factor(id) ) )
#GOexp <- subset( Dgo, complete.cases(rt) & cond == "exp" ) %>% mutate( id = as.integer( as.factor(id) ) ) 
#SRcon <- subset( Dsr, complete.cases(rt) & cond == "ctrl" & (rt>ssd) ) %>% mutate( id = as.integer( as.factor(id) ) )
#SRexp <- subset( Dsr, complete.cases(rt) & cond == "exp" & (rt>ssd) ) %>% mutate( id = as.integer( as.factor(id) ) )
#NAcon <- subset( Dna, cond == "ctrl" ) %>% mutate( id = as.integer( as.factor(id) ) )
#NAexp <- subset( Dna, cond == "exp" ) %>% mutate( id = as.integer( as.factor(id) ) )


# DCM MODEL ----

# read the files
source( here("_dmc","DMC_190819","dmc","dmc.R") )
load_model("EXG-SS","exgSS.R")
