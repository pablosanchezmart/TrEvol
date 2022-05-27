#' Plot variance-covariance networks
#'
#' @param covariance.results (data frame) correlations output of correlationsTraits$covarianceResults
#' @param variance.results (data frame) phylogenetic signal as reported by phylogeneticSignalTraits$varianceResults
#' @param variance.type (character) variable name of the phylogenetic signal as reported in varianceResults
#' @param covariance.type (character) variable name of the correlation as reported iun correlations
#' @param group.variables (data frame) character vector as long as the number of variables indicating the variable name in the first column and the group of the variable in the second.
#' @param order.variables (character) order of the variables to plot
#' @param edge.label (logical) if TRUE, edge labels show correlation coefficients
#' @param layout "spring" or "circular"
#' @param only.significant (logical) if TRUE, only significant correlations are plotted
#' @param not.show.variables (character) vector of variables present in correlations data but exclued from the plot
#' @param threshold (numeric) correlations to be shown (e.g., those with a absolute value higher than 0.1)
#' @param label.size (numeric) size of the node labels
#' @param node.label (data frame) character vector as long as the number of variables indicating the variable name in the first column and the new name, which will be shown in the plot.
#'
#' @return
#' @export
#'
#' @examples
plotNetwork <- function (covariance.results, variance.results = NULL, variance.type = NULL,
                         covariance.type, group.variables = NULL, order.variables = NULL, edge.label = F,
                         layout = "spring", only.significant = T, not.show.variables = NULL,
                         threshold = 0, label.size = 0.8, node.label = NULL, plot.metrics = T)
{
  # set non significant values to zero if needed
  if (only.significant) {
    covariance.results[which(covariance.results[, paste0("Pvalue_", covariance.type)] > 0.05), covariance.type] <- 0
  }
  correlations <- covariance.results[, c("Trait_1", "Trait_2", covariance.type)]
  vars <- unique(c(correlations$Trait_1, correlations$Trait_1))
  correlation.matrix <- matrix(ncol = length(vars), nrow = length(vars),
                               0)
  colnames(correlation.matrix) <- vars
  rownames(correlation.matrix) <- vars
  for (i in 1:length(correlations$Variable1)) {
    var1 <- correlations$Variable1[i]
    var2 <- correlations$Variable2[i]
    correlation.matrix[var1, var2] <- correlations[i, covariance.type]
    correlation.matrix[var2, var1] <- correlations[i, covariance.type]
    correlation.matrix[var1, var1] <- 0
    correlation.matrix[var2, var2] <- 0
  }
  if (!is.null(not.show.variables) && any(not.show.variables %in% vars)) {
    vars <- vars[-which(vars %in% not.show.variables)]
    correlation.matrix <- correlation.matrix[vars, vars]
  }
  if (!is.null(group.variables)) {
    if (all(vars %in% group.variables[, 1])) {
      group.variables <- group.variables[which(group.variables[, 1] %in% vars),
      ]
      correlation.matrix <- correlation.matrix[group.variables[, 1], group.variables[, 1]]
      vars <- group.variables[, 1]
    } else {
      stop("Not all variables of interest are included in group.variables argument.")
    }
  }
  if (!is.null(order.variables)) {
    if (all(vars %in% order.variables)) {
      vars <- order.variables[which(order.variables %in% vars)]
      correlation.matrix <- correlation.matrix[vars, vars]
    } else {
      stop("Not all variables of interest are included in order.variables argument.")
    }
  }
  corGraph <- igraph::graph.adjacency(abs(as.matrix(correlation.matrix)),
                                      mode = "lower", weighted = T)
  diameter <- round(igraph::diameter(corGraph, directed = F),
                    3)
  beta_connectivity <- round(igraph::ecount(corGraph)/igraph::vcount(corGraph),
                             3)
  transitivity <- round(igraph::transitivity(corGraph, isolates = "zero",
                                             type = "global"), 3)
  networkMetrics <- data.frame(correlation_type = covariance.type,
                               diameter = diameter, beta_connectivity = beta_connectivity,
                               transtiivity = transitivity)
  rownames(node.label) <- node.label[, 1]
  node.label <- node.label[vars, 2]
  if (!is.null(variance.results)) {
    rownames(variance.results) <- variance.results[, 1]
    if (!is.null(variance.results) && all(vars %in% rownames(variance.results))) {
      ps.vars <- variance.results[vars, variance.type]
    } else {
      stop("Not all variables of interst are present in the phylogenetic signal dataframe")
    }
  } else {
    ps.vars <- rep(NULL, length(correlation.matrix[, 1]))
  }
  p <- qgraph::qgraph(correlation.matrix, layout = layout,
                      vsize = 5, vsize2 = 1, esize = 10 * max(correlation.matrix),
                      palette = "pastel", negDashed = T, borders = T, legend = F,
                      vTrans = 180, fade = F, aspect = T, legend.cex = 0.25,
                      edge.labels = edge.label, edge.label.cex = 2 * length(vars)/length(vars),
                      labels = node.label, label.cex = label.size, label.scale = F,
                      node.label.offset = c(0.5, -2), pie = ps.vars, pieBorder = 1,
                      DoNotPlot = F, threshold = threshold, groups = group.variables[, 2])

  if(plot.metrics){
    graphics::text(x = 0.5, y = -1, labels = paste0("Connectivity = ",
                                                    round(networkMetrics$beta_connectivity, 2), "\n", "Transitivity = ",
                                                    round(networkMetrics$transtiivity, 2)), adj = 0, cex = 0.6)
  }
  return(p)
}
