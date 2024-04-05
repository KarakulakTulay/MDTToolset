#' Title
#'
#' @param protein ENSP id of protein of interest
#' @param iso_network isoform interaction network
#' @return plot
#' @import igraph
#' @importFrom RColorBrewer brewer.pal
#' @export

visualize_network <- function(protein, iso_network) {
  # create data:
  links <- data.frame(
    source=c(protein),
    target=c(regmatches(iso_network[iso_network$ENSP == protein, 'MissInts'] , gregexpr("ENSP\\d+", iso_network[iso_network$ENSP == protein, 'MissInts'] ))[[1]], regmatches(iso_network[iso_network$ENSP == protein, 'ExistInts'] , gregexpr("ENSP\\d+", iso_network[iso_network$ENSP == protein, 'ExistInts'] ))[[1]]),
    importance=(sample(1, 1, replace=T))
  )
  nodes <- data.frame(
    name=c(protein, regmatches(iso_network[iso_network$ENSP == protein, 'MissInts'] , gregexpr("ENSP\\d+", iso_network[iso_network$ENSP == protein, 'MissInts'] ))[[1]], regmatches(iso_network[iso_network$ENSP == protein, 'ExistInts'] , gregexpr("ENSP\\d+", iso_network[iso_network$ENSP == protein, 'ExistInts'] ))[[1]]),
    carac=c('Main', rep('Int Loss', length(regmatches(iso_network[iso_network$ENSP == protein, 'MissInts'] , gregexpr("ENSP\\d+", iso_network[iso_network$ENSP == protein, 'MissInts'] ))[[1]])), rep('Int Kept', length(regmatches(iso_network[iso_network$ENSP == protein, 'ExistInts'] , gregexpr("ENSP\\d+", iso_network[iso_network$ENSP == protein, 'ExistInts'] ))[[1]])))
  )

  # Turn it into igraph object
  network <- igraph::graph_from_data_frame(d=links, vertices=nodes, directed=F)

  # Make a palette of 3 colors
  coul  <- RColorBrewer::brewer.pal(3, "Dark2")

  # Create a vector of color
  my_color <- coul[as.numeric(as.factor(V(network)$carac))]

  # Make the plot
  network_plot <- plot(network, vertex.color=my_color, vertex.size=20, vertex.label.dist = 1 , edge.width = 4, vertex.label.cex = 0.75, vertex.label.degree=-pi/2)

  # Add a legend
  legend("bottomleft", legend=levels(as.factor(V(network)$carac))  , col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 0.75, text.col=coul , horiz = FALSE, inset = c(0.01, 0.01))

  return(network_plot)
}
