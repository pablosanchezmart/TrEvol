#' Plot variance-covariance networks
#'
#' @param correlations (data frame) correlations output of correlationsTraits$correlations.rslts
#' @param phylogenetic.signal (data frame() phylogenetic signal as reported by phylogeneticSignalTraits$phylogenetic.signal.rslts
#' @param phyloSignal.name (character) variable name of the phylogenetic signal as reported in phylogenetic.signal
#' @param correlation.type (character) variable name of the correlation as reported iun correlations
#' @param gr_vars (data frame) character vector as long as the number of variables indicating the variable name in the first column and the group of the variable in the second.
#' @param order_vars (character) order of the variables to plot
#' @param edgeLab (logical) if TRUE, edge labels show correlation coefficients
#' @param layout "spring" or "circular"
#' @param onlySignificant (logical) if TRUE, only significant correlations are plotted
#' @param notShowCors (character) vector of variables present in correlations data but exclued from the plot
#' @param threshold (numeric) correlations to be shown (e.g., those with a absolute value higher than 0.1)
#' @param labelSize (numeric) size of the node labels
#' @param nodeLab (data frame) character vector as long as the number of variables indicating the variable name in the first column and the new name, which will be shown in the plot.
#'
#' @return
#' @export
#'
#' @examples
plotNetwork <- function(correlations,
                                     phylogenetic.signal = NULL,
                                     phyloSignal.name = NULL,
                                     correlation.type,
                                     gr_vars = NULL,
                                     order_vars = NULL,
                                     edgeLab = F,
                                     layout = "spring",
                                     onlySignificant = T,
                                     notShowCors = NA,
                                     threshold = 0,
                                     labelSize = 0.8,
                                     nodeLab) {

  # Correlation matrix
  if(onlySignificant){
    correlations[which(correlations[ , paste0(correlation.type, "_pvalue")] > 0.05), correlation.type] <- 0 # non significant correlation coefficients set to zero
  }

  correlations <- correlations[, c("Variable1", "Variable2", correlation.type)]

  vars <- unique(c(correlations$Variable1, correlations$Variable2))

  correlation.matrix <- matrix(ncol = length(vars), nrow = length(vars), 0)
  colnames(correlation.matrix) <- vars
  rownames(correlation.matrix) <- vars

  for(i in 1:length(correlations$Variable1)){
    var1 <- correlations$Variable1[i]
    var2 <- correlations$Variable2[i]
    correlation.matrix[var1, var2] <- correlations[i, correlation.type]
    correlation.matrix[var2, var1] <- correlations[i, correlation.type]
    correlation.matrix[var1, var1] <- 0
    correlation.matrix[var2, var2] <- 0
  }

  if(!is.na(notShowCors) && any(notShowCors %in% vars)){
    vars <- vars[-which(vars %in% notShowCors)]
    correlation.matrix <- correlation.matrix[vars, vars]
  }

  if(!is.null(gr_vars)){
    gr_vars <- gr_vars %>% filter(variable %in% vars)
    correlation.matrix <- correlation.matrix[gr_vars[, 1], gr_vars[, 1]]
    vars <- gr_vars[, 1]
  }

  if(!is.null(order_vars)){
    vars <- order_vars[which(order_vars %in% vars)]
    correlation.matrix <- correlation.matrix[vars, vars]
  }

  # Network metrics

  corGraph <- igraph::graph.adjacency(abs(as.matrix(correlation.matrix)), mode = "lower", weighted = T)
  # Network metrics calculation
  diameter <- round(diameter(corGraph, directed = F), 3)
  beta_connectivity <- round(igraph::ecount(corGraph)/igraph::vcount(corGraph), 3)
  transitivity <- round(igraph::transitivity(corGraph, isolates = "zero", type = "global"), 3) # type = "weighted"

  networkMetrics <- data.frame(correlation_type = correlation.type,
                               "diameter" = diameter,
                               "beta_connectivity" = beta_connectivity,
                               "transtiivity" = transitivity)

  ## Node labels
  rownames(nodeLab) <- nodeLab$variable
  nodeLab <- nodeLab[vars, "name"]

  ## Phylogenetic signal

  rownames(phylogenetic.signal) <- phylogenetic.signal[, 1]
  ps.vars <- phylogenetic.signal[vars, phyloSignal.name]

  p <- qgraph::qgraph(correlation.matrix, layout = layout,
              vsize = 5, vsize2 = 1, esize = 10*max(correlation.matrix),
              palette = "pastel",  negDashed = T, borders =  T, legend = F, vTrans = 180, fade = T, aspect = T, legend.cex = 0.25,# posCol = "#377eb8", negCol = "#e41a1c",
              edge.labels = edgeLab, edge.label.cex = 2*length(nodeLab)/length(nodeLab),
              labels = nodeLab, label.cex = labelSize, label.scale = F,
              node.label.offset = c(0.5, -2),
              pie = ps.vars, pieBorder = 1, # pieColor = "black",
              DoNotPlot = F,
              threshold = threshold,
              groups = gr_vars[, 2]
  )
  graphics::text(x=0.5, y=-1, labels= paste0("Connectivity = ", round(networkMetrics$beta_connectivity, 2), "\n",
                                   "Transitivity = ", round(networkMetrics$transtiivity, 2)), adj = 0,
       cex = 0.6)

  return(p)
}
