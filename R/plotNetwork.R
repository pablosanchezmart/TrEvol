#' Plot variance-covariance networks
#'
#' @param covariance.results (data frame) correlations output of computeCovariancePartition()$covarianceResults
#' @param variance.results (data frame) phylogenetic signal as reported by computeVariancePartition()$varianceResults
#' @param variance.type (character) variable name of the phylogenetic signal as reported in varianceResults
#' @param covariance.type (character) variable name of the correlation as reported in correlations
#' @param group.variables (data frame) character vector as long as the number of variables indicating the variable name in the first column and the group of the variable in the second.
#' @param order.variables (character) order of the variables to plot
#' @param edge.label (logical) if TRUE, edge labels show correlation coefficients
#' @param layout "spring" or "circular"
#' @param only.significant (logical) if TRUE, only significant correlations are plotted
#' @param not.show.variables (character) vector of variables present in correlations data but exclued from the plot
#' @param threshold (numeric) correlations to be shown (e.g., those with a absolute value higher than 0.1)
#' @param label.size (numeric) size of the node labels
#' @param node.label (data frame) character vector as long as the number of variables indicating the variable name in the first column and the new name, which will be shown in the plot.
#' @param networkMetrixTextSize (numeric) Size of the text reporting network metrics
#' @param displayDegreeAsNodeSize (logical) If true, node size is proportional to its degree.
#'
#' @return
#' @export
#'
#' @examples
plotNetwork <- function (covariance.results, variance.results = NULL, variance.type = NULL,
                         covariance.type, group.variables = NULL, order.variables = NULL, edge.label = F,
                         layout = "spring", only.significant = T, not.show.variables = NULL,
                         threshold = 0, label.size = 0.8, node.label = NULL, plot.metrics = T,
                         networkMetrixTextSize = 1, displayDegreeAsNodeSize = T)
{

  ### Data preparation ####

  # set non significant values to zero if needed
  if (only.significant) {
    covariance.results[which(covariance.results[, paste0("Pvalue_", covariance.type)] > 0.05), covariance.type] <- 0
  }

  if (threshold > 0) {
    covariance.results[which(abs(covariance.results[, covariance.type]) < threshold), covariance.type] <- 0
  }
  correlations <- covariance.results[, c("Trait_1", "Trait_2", covariance.type)]
  vars <- unique(c(correlations$Trait_1, correlations$Trait_2))
  correlation.matrix <- matrix(ncol = length(vars), nrow = length(vars), 0)
  colnames(correlation.matrix) <- vars
  rownames(correlation.matrix) <- vars
  for (i in 1:length(correlations$Trait_1)) {
    var1 <- correlations$Trait_1[i]
    var2 <- correlations$Trait_2[i]
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

  ### Network metrics ####

  ## Graph
  corGraph <- igraph::graph.adjacency(abs(as.matrix(correlation.matrix)),
                                      mode = "lower", weighted = T)

  ## Edge density (ED)
  # ED describes the density of the connected edges between nodes in a network, that is, the proportion of actual connections
  # among traits out of all possible connections.

  EdgeDensity <- (2* igraph::ecount(corGraph)) / (igraph::vcount(corGraph) * (igraph::vcount(corGraph) -1))

  ## Diameter (D) and average path length (AL)
  # D is the maximum shortest distance between any two connected node traits in the network, and AL is the mean shortest path between all node traits
  # in the network. PTNs with higher D and AL have greater overall independence among traits.

  diameter <- round(igraph::diameter(corGraph, directed = F), 3)
  averagePathLength <- igraph::mean_distance(corGraph)

  ## Average clustering (AC)
  # AC is the average of the clustering coefficients of all traits in PTNs. PTNs with higher AC are more extensively divided into
  # several different components.

  averageClusteringCoefficient <- mean(round(igraph::transitivity(corGraph, isolates = "zero", type = "local"), 3))

  ## Node degree (k)

  # k is the number of edges that connect a focal node trait to other nodes. Plant traits that have a high k can be considered as
  # overall hub traits.

  degree <- igraph::degree(corGraph)


  networkMetrics <- data.frame(correlation_type = covariance.type,
                               EdgeDensity = EdgeDensity,
                               diameter = diameter,
                               averagePathLength = averagePathLength,
                               averageClusteringCoefficient = averageClusteringCoefficient)

  rownames(node.label) <- node.label[, 1]
  node.label <- node.label[vars, 2]

  ### Variance results (node pie chart) ####

  if (!is.null(variance.results)) {
    rownames(variance.results) <- variance.results[, 1]
    if (!is.null(variance.results) && all(vars %in% rownames(variance.results))) {
      ps.vars <- variance.results[vars, variance.type]
      ps.vars[which(ps.vars < 0)] <- 0
      ps.vars[which(ps.vars > 1)] <- 1
    } else {
      stop("Not all variables of interest are present in the phylogenetic signal dataframe")
    }
  } else {
    ps.vars <- rep(NULL, length(correlation.matrix[, 1]))
  }

  ### Plot ####

  if(displayDegreeAsNodeSize){
    nodeSize <- 5 + degree * 2
    nodeSize <- ifelse(nodeSize > 10, 10, nodeSize)
    nodeSize <- ifelse(nodeSize < 5, 5, nodeSize)
  } else {
    nodeSize <- 8
  }

  p <- qgraph::qgraph(correlation.matrix, layout = layout,
                      vsize = nodeSize, vsize2 = nodeSize, esize = 10 * max(abs(correlation.matrix)),
                      palette = "pastel", negDashed = T, borders = T, legend = F,
                      vTrans = 180, fade = F, aspect = T, legend.cex = 0.25,
                      edge.labels = edge.label, edge.label.cex = 2 * length(vars)/length(vars),
                      labels = node.label, label.cex = label.size, label.scale = F,
                      node.label.offset = c(0.5, -3), pie = ps.vars, pieBorder = 1,
                      DoNotPlot = F, groups = group.variables[, 2], threshold = threshold)

  if(plot.metrics){
    graphics::text(x = 0.5, y = -1, labels = paste0("ED = ", round(networkMetrics$EdgeDensity, 2), "\n",
                                                    "D = ", round(networkMetrics$diameter, 2), "\n",
                                                    "AL = ", round(networkMetrics$averagePathLength, 2), "\n",
                                                    "AC = ", round(networkMetrics$averageClusteringCoefficient, 2)),
                   adj = 0, cex = networkMetrixTextSize)
  }
  return(p)
}
