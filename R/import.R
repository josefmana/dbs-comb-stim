#
# In this script I prepare in-house functions used for data import
#

#
# COUNTERBALANCING ----
generate_counterbalancing <- function(.file) {
  
  # read the file and add column with counterbalancing (round with added medium-frequency stimulation)
  data <- .file %>% mutate(LFS = if_else(is.na(round_1), 2, 1) )
  cb0 <- lapply(1:2, function(i) subset(data, LFS == i)$id) # list patients in each condition
  
  # continue with each subset of participants separately
  # compute all possibilities where half of patients change condition
  # one of which are from LFS == 2, and three of which are from LFS == 1
  cb_change0 <- lapply(
    
    cb0, function(id)
      
      lapply( set_names( data$id[data$id %in% id] ), function(i) 0:1) %>% # list of patients' ids
      do.call(expand_grid, . ) %>% # get all combinations
      filter( rowSums(.) == 2) %>% # keep only those that sum to 2 (i.e., 2/3 in LFS == 2, 2/5 in LFS == 1)
      mutate_all(~ 1-.x) # flip change signs (leading to 1/3 changing from LFS == 2, 3/5 changing from LFS == 1)
    
  )
  
  # combine possibilities of valid counterbalancing change fro LFS == 1 and LFS == 2 groups
  possibilities <- expand_grid( `1` = 1:nrow(cb_change0[[1]]), `2` = 1:1:nrow(cb_change0[[2]]) )
  
  # list possible combinations of 'stay' (0) vs 'change' (1) counterbalancing
  cb_change <- sapply(
    
    1:nrow(possibilities),
    function(i) cbind.data.frame(
      
      cb_change0[[1]][as.numeric(possibilities[i, 1]), ], # originally LFS == 1 patients
      cb_change0[[2]][as.numeric(possibilities[i, 2]), ]  # originally LFS == 2 patients
      
    ) %>% unlist() # unlist to allow for colSums later
  )
  
  # extract counterbalancing options for the re-test
  cb_new <-
    cb_change %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      across(
        everything(),
        ~ if_else(
          .x == 0,
          which( sapply(1:2, function(i) cur_column() %in% cb0[[i]] ) ),
          which( !sapply(1:2, function(i) cur_column() %in% cb0[[i]] ) )
        )
      )
    ) %>%
    t()
    
  # some checks
  cb_change %>% colSums() # change, should be 4 (i.e., half of the participants change counterbalancing)
  bcb_new %>% as_tibble() %>% mutate_all(~ .x-1) %>% colSums() %>% as.numeric() # LFS == 2, should be 5 in the re-test
  cb_new[cb0[[1]] ,] %>% as_tibble() %>% mutate_all(~ .x-1) %>% colSums() %>% as.numeric() # number of LFS == 1 -> LFS == 2 changes, should be 3
  cb_new[cb0[[2]] ,] %>% as_tibble() %>% mutate_all(~ .x-1) %>% colSums() %>% as.numeric() # number of LFS == 2 -> LFS == 2 'changes', should be 2
  
}
