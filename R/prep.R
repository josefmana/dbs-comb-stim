#
# In this script I prepare in-house functions used for data import
#

#
# COUNTERBALANCING ----
generate_counterbalancing <- function(.file) {
  
  data <- .file %>% mutate(LFS = if_else(is.na(round_1), 2, 1) ) # read the file and add column with counterbalancing (round with added medium-frequency stimulation
  cb0 <- lapply(1:2, function(i) subset(data, LFS == i)$id) # list patients in each condition
  
  # continue with each subset of participants separately
  # compute all possibilities wherein half of patients change condition
  # one of which is from LFS == 2 to LFS == 1,
  # and three of which are from LFS == 1 to LFS == 2
  cb_change_list <- lapply(
    
    cb0, function(id)
      
      lapply( set_names( data$id[data$id %in% id] ), function(i) 0:1) %>% # list of patients' ids
      do.call(expand_grid, . ) %>% # get all combinations
      filter( rowSums(.) == 2) %>% # keep only those that sum to 2 (i.e., 2/3 in LFS == 2, 2/5 in LFS == 1)
      mutate_all(~ 1-.x) # flip change signs (leading to 1/3 changing from LFS == 2, 3/5 changing from LFS == 1)
    
  )
  
  # combine possibilities (row combinations) of valid counterbalancing change fro LFS == 1 and LFS == 2 groups
  possibilities <- expand_grid( `1` = 1:nrow(cb_change_list[[1]]), `2` = 1:1:nrow(cb_change_list[[2]]) )
  
  # list possible combinations of 'stay' (0) vs 'change' (1) counterbalancing
  cb_change <- sapply(
    
    1:nrow(possibilities),
    function(i) cbind.data.frame(
      
      cb_change_list[[1]][as.numeric(possibilities[i, 1]), ], # originally LFS == 1 patients
      cb_change_list[[2]][as.numeric(possibilities[i, 2]), ]  # originally LFS == 2 patients
      
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
          condition = .x == 0,
          true = which( sapply(1:2, function(i) cur_column() %in% cb0[[i]] ) ),
          false = which( !sapply(1:2, function(i) cur_column() %in% cb0[[i]] ) )
        )
      )
    ) %>%
    t()
  
  # some checks
  ctrl <- list(
    
    changed = list(
      observed = cb_change %>% colSums(),
      desired = "Exacly 4 participants ought to change their counterbalancing in. the re-test."
    ),
    lfs2 = list(
      observed = cb_new %>% as_tibble() %>% mutate_all(~ .x-1) %>% colSums() %>% as.numeric(),
      desired = "Exactly 5 particinats ought to be in the LFS == 2 condition in the re-test."
    ),
    lfs1_to_lfs2 = list(
      observed = cb_new[cb0[[1]] ,] %>% as_tibble() %>% mutate_all(~ .x-1) %>% colSums() %>% as.numeric(),
      desired = "Exactly 3 participants ought to change from LFS == 1 in the original to LFS == 2 in the re-test."
    ),
    lfs2_to_lfs2 = list(
      observed = cb_new[cb0[[2]] ,] %>% as_tibble() %>% mutate_all(~ .x-1) %>% colSums() %>% as.numeric(),
      desired = "Exactly 2 participants ought to not change from LFS == 2 in the original to the re-test."
    )
    
  )
  
  # return resulting counterbalancing options and their checks
  return(
    list(
      options = cb_new,
      checks = ctrl
    )
  )
  
}
