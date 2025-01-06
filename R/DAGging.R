#
# This script contains functions for causal assumption representation.
#

make_dag <- function(plot = T, save = T) {
  
  # prepare names
  nms <- data.frame(
    
    label = c("Frequency", "Individual", "Contact\nLocation", "Motor\nSpeed", "Processing\nSpeed", "Response\nInhibition", "Finger\nTapping", "Response\nTime", "Trial\nType"),
    name = c("Freq", "Indi", "CLoc", "MS", "PS", "RI", "FT", "RT", "TT"),
    x = rep(1:3, 3),
    y = c( 3, 4, 3, rep(2, 3), rep(1, 3) )
    
  ) %>% mutate(
    
    latent = if_else(name %in% c("MS","PS","RI") , T, F),
    colour = if_else(latent == T, "black", "white")
    
  )
  
  # list latent variables separately
  lvs <- with(nms, name[latent == T])
  
  # prepare the DAG
  dag <- dagify(
    
    RT ~ MS + PS + RI + TT,
    FT ~ MS,
    MS ~ Freq + CLoc + Indi,
    PS ~ Freq + CLoc + Indi,
    RI ~ Freq + CLoc + Indi,
    CLoc ~ Indi,
    
    latent = lvs,
    coords = nms[ , c("name","x","y")]
    
  ) %>%
    
    tidy_dagitty() %>%
    arrange(name) %>%
    mutate( latent = if_else(name %in% lvs, "1", "0") )
  
  # plot the DAG
  dagplt <- dag %>%
    
    ggplot() +
    aes(x = x, y = y, xend = xend, yend = yend, colour = latent) +
    
    geom_dag_point(size = 20, fill = "white", shape = 21, stroke = 1.1) +
    geom_dag_edges(arrow_directed = grid::arrow(length = grid::unit(13, "pt"), type = "open"), edge_width = .89) +
    scale_colour_manual( values = c("white", "black") ) +
    geom_dag_text(
      label = arrange(nms, name)$name,
      color = "black",
      size = 6
    ) +
    
    theme_dag() +
    theme(legend.position = "none")
  
  # save & return it
  if(save == T) ggsave(plot = dagplt, filename = here("DAG.jpg"), dpi = 300, width = 7, height = 7)
  if(plot == T) return(dagplt) else return(dag)
  

}

