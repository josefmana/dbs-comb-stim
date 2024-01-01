# This is a script that takes reports from 'reps' folder to their respective folders labelled according to output format.

rm( list = ls() ) # clean environment
library(here) # load packages needed

# list formats of interest
f <- c("html","pdf","docx")

# prepare folders for html, docx and pdf outputs
sapply( f, function(i) if( !dir.exists(i) ) dir.create(i) )

# look through 'reps' folder and move reports
for( i in f ) {
  
  # list all reports in format i
  files <- list.files("reps")[ grepl( paste0(".",i), list.files("reps") ) ]
  
  for ( j in files ) {
    
    # move it
    file.copy( from = here("reps",j), to = here(i,j), overwrite = T )
    file.remove( here("reps",j) )
    
  }
  
}
