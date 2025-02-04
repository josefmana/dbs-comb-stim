#
# This script contains functions for causal assumption representation.
#

#
# DRAW A DAG ----
make_dag <- function(plot = T, save = T) {
  
  # prepare names
  nms <- data.frame(
    
    label = c(
      "Frequency",
      "Individual",
      "Contact\nLocation",
      "Speech\nProduction",
      "Processing\nSpeed",
      "Response\nInhibition",
      "Affect",
      "Words",
      "Response\nTime",
      "STAI",
      "Probe[VF]",
      "Mode[VF]",
      "Trial\nType[SSRT]",
      "Version[STAI]",
      "Item[STAI]"
    ),
    name  = c(
      "Freq",
      "Indi",
      "CLoc",
      "SP",
      "PS",
      "RI",
      "Aff",
      "VF",
      "RT",
      "STAI",
      "Probe",
      "Mode",
      "TT",
      "X",
      "Item"
    ),
    x = c(
      2:4,                                   # Individual & experimental manipulations
      seq(from = 1, to = 5, length.out = 4), # Latent variables
      2:4,                                   # Observed variables
      1:5                                    # Observed variables' structure
    ),
    y = c(
      4, 5, 4,   # Individual & experimental manipulations
      rep(3, 4), # Latent variables
      rep(2, 3), # Observed variables
      rep(1, 5)  # Observed variables' structure
    )
    
  ) %>% mutate(
    
    latent = if_else(name %in% c("SP", "PS", "RI", "Aff") , T, F),
    colour = if_else(latent == T, "black", "white")
    
  )
  
  # list latent variables separately
  lvs <- with(nms, name[latent == T])
  
  # prepare the DAG
  dag <- dagify(
    
    CLoc ~ Indi,
    
    SP  ~ Freq + CLoc + Indi,
    PS  ~ Freq + CLoc + Indi,
    RI  ~ Freq + CLoc + Indi,
    Aff ~ Freq + CLoc + Indi,
    
    VF   ~ PS + SP + Mode + Probe,
    RT   ~ PS + RI + TT,
    STAI ~ Aff + X + Item + Indi,
    
    Probe ~ Mode,
    Item  ~ X,
    
    latent = lvs,
    coords = nms[ , c("name","x","y")]
    
  ) %>%
    
    tidy_dagitty() %>%
    arrange(name) %>%
    mutate(
      latent = if_else(
        condition = name %in% lvs,
        true      = "1",
        false     = "0"
      ),
      curve  = if_else(
        condition = is.na(direction),
        true      = NA,
        false     = case_when(
          name == "Indi" & to == "STAI" ~ -0.24,
          name == "Indi" & to == "SP"   ~ -0.33,
          name == "Indi" & to == "Aff"  ~  0.33,
          .default = 0
        )
      )
    )
  
  # plot the DAG
  dagplt <- dag %>%
    
    ggplot() +
    aes(x = x, y = y, xend = xend, yend = yend, colour = latent) +
    
    geom_dag_point(size = 20, fill = "white", shape = 21, stroke = 1.1) +
    geom_dag_edges_arc(
      curvature  = na.omit(dag$data$curve),
      arrow      = grid::arrow(length = grid::unit(13, "pt"), type = "open"),
      edge_width = .89
    ) +
    scale_colour_manual( values = c("white", "black") ) +
    geom_dag_text(
      label = arrange(nms, name)$name,
      color = "black",
      size = 6
    ) +
    
    theme_dag() +
    theme(legend.position = "none")
  
  # save & return it
  if(save == T) ggsave(plot = dagplt, filename = here("DAG.jpg"), dpi = 300, width = 9, height = 11)
  if(plot == T) return(dagplt) else return(dag)

}

