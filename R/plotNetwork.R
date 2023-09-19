#' Plot variance-covariance networks
#'
#' Plot variance-covariance partition obtained using the function computeVarianceCovariancePartition of this package.
#'
#' @param variance_results (*data frame*). Variance partition as reported by computeVarianceCovariance(). It can be found in "results_object$varianceResults.
#' @param correlation_results (*data frame*). Correlations output of computeVarianceCovariance(). It can be found in "results_object$covarianceResults".
#' @param variance_type (*character*). Name of the portion of variance to report in the nodes of the plot. When environmental variable is not included in computeVarianceCovariance(), it can be set as "phylogenetic_variance" or "non_phylogenetic_variance".
#' If environmental variable is included, this argument can be set as "non_attributed_phylogenetic_variance", "environmental_phylogenetic_variance", "labile_environmental_variance" and "residual_variance".
#' @param correlation_type (*character*). Name of the portion of the correlation to report in the edges of the plot. When environmental variable is not included in computeVarianceCovariance(), it can be set as "total_correlation", "phylogenetic_correlation" or "non_phylogenetic_correlation".
#' If environmental variable is included, this argument can be set as "non_attributed_phylogenetic_correlation", "environmental_phylogenetic_correlation", "labile_environmental_correlation" and "residual_correlation".
#' @param group_variables (*data frame*). Data frame wich a first column indicating the variables to be plotted included in correlation_results and varaince_results and the group of the variable in the second column. Variables of the same group will be ploted with the same colour.
#' @param order_variables (*character*). Character indicating the variables names to be plotted in the desired order. The variables will be plotted in the desired order clockwise, when layout is set to be circular.
#' @param edge_label (*logical*). If TRUE, edge labels show correlation coefficients. Default set to FALSE.
#' @param layout (*"spring"*, or *"circular"*). Network layout (default as circular). See qgraph documentation for further information.
#' @param only_significant (*logical*). if TRUE (default), only significant correlations are represented.
#' @param exclude_variables (*character*). Vector of variables present in correlations data but excluded from the plot. It can be useful to plot subsets of the results.
#' @param correlation_threshold (*numeric*). Only correlation with an absolute value higher than the threshold will be plotted. When set to zero, there is no threshold (default).
#' @param label_size (*numeric*). Control the size of the node labels.
#' @param node_label (*data frame*). Data frame indicating the variable name in the first column and the new name, which will be shown in the plot, in the second column.
#' @param plot_metrics (*logical*). If true, network metrics are displayed.
#' @param network_metrics_text_size (*numeric*). Size of the text reporting network metrics.
#' @param degree_as_node_size (*logical*) If TRUE (default), node size is proportional to the degree of the variable (number of edges connected to a given node representing a variable).
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate example data
#' simulated_traits.data <- simulateDataSet()
#'
#' # Compute variance-covariance structure for simulated traits using default parameters
#' variance_correlation_results <- computeVarianceCovariancePartition(
#' traits = c("phylo_G1_trait1", "phylo_G1_trait2"),
#' environmental_variable = "phylo_G1_env",
#' phylogeny = simulated_traits.data$phylogeny
#' )
#'
#' # plot network
#'
#' plotNetwork(
#' correlation_results = variance_correlation_results$covarianceResults,
#' variance_results = variance_correlation_results$varianceResults,
#' variance_type = environmental_phylogenetic_variance,
#' correlation_type = environmental_phylogenetic_covariance
#' )
#' }
plotNetwork <- function (variance_results = NULL,
                         correlation_results = NULL,
                         variance_type = NULL,
                         correlation_type = NULL,
                         group_variables = NULL,
                         order_variables = NULL,
                         edge_label = F,
                         layout = "circular",
                         only_significant = T,
                         exclude_variables = NULL,
                         correlation_threshold = 0,
                         label_size = 0.8,
                         node_label = NULL,
                         plot_metrics = T,
                         network_metrics_text_size = 1,
                         degree_as_node_size = T)
{

  # Arguments
  if(is.null(variance_results)){
    stop("Specify variance_results argument")
  }

  if(is.null(correlation_results)){
    stop("Specify correlation_results argument")
  }

  if(is.null(variance_type)){
    stop("Specify variance_type argument")
  }

  if(is.null(correlation_type)){
    stop("Specify correlation_type argument")
  }

  if(is.null(correlation_type)){
    stop("Specify correlation_type argument")
  }

  ### Data preparation ####

  # set non significant values to zero if needed
  if (only_significant) {
    correlation_results[which(correlation_results[, paste0("p_value_", correlation_type)] > 0.05), correlation_type] <- 0
  }

  if (correlation_threshold > 0) {
    correlation_results[which(abs(correlation_results[, correlation_type]) < correlation_threshold), correlation_type] <- 0
  }
  correlations <- correlation_results[, c("trait_1", "trait_2", correlation_type)]
  vars <- unique(c(correlations$trait_1, correlations$trait_2))
  correlation.matrix <- matrix(ncol = length(vars), nrow = length(vars), 0)
  colnames(correlation.matrix) <- vars
  rownames(correlation.matrix) <- vars
  for (i in 1:length(correlations$trait_1)) {
    var1 <- correlations$trait_1[i]
    var2 <- correlations$trait_2[i]
    correlation.matrix[var1, var2] <- correlations[i, correlation_type]
    correlation.matrix[var2, var1] <- correlations[i, correlation_type]
    correlation.matrix[var1, var1] <- 0
    correlation.matrix[var2, var2] <- 0
  }
  if (!is.null(exclude_variables) && any(exclude_variables %in% vars)) {
    vars <- vars[-which(vars %in% exclude_variables)]
    correlation.matrix <- correlation.matrix[vars, vars]
  }
  if (!is.null(group_variables)) {
    if (all(vars %in% group_variables[, 1])) {
      group_variables <- group_variables[which(group_variables[, 1] %in% vars),
      ]
      correlation.matrix <- correlation.matrix[group_variables[, 1], group_variables[, 1]]
      vars <- group_variables[, 1]
    } else {
      stop("Not all variables of interest are included in group_variables argument.")
    }
  }
  if (!is.null(order_variables)) {
    if (all(vars %in% order_variables)) {
      vars <- order_variables[which(order_variables %in% vars)]
      correlation.matrix <- correlation.matrix[vars, vars]
    } else {
      stop("Not all variables of interest are included in order_variables argument.")
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

  # diameter <- round(igraph::diameter(corGraph, directed = F), 3)
  diameter <- max(abs(correlation.matrix))

  mean_path_length_matrix <- correlation.matrix
  mean_path_length_matrix[mean_path_length_matrix == 0] <- NA
  # averagePathLength <- igraph::mean_distance(corGraph)
  averagePathLength <- mean(abs(mean_path_length_matrix), na.rm = T)

  ## Average clustering (AC)
  # AC is the average of the clustering coefficients of all traits in PTNs. PTNs with higher AC are more extensively divided into
  # several different components.

  averageClusteringCoefficient <- mean(round(igraph::transitivity(corGraph, isolates = "zero", type = "local"), 3))

  ## Node degree (k)

  # k is the number of edges that connect a focal node trait to other nodes. Plant traits that have a high k can be considered as
  # overall hub traits.

  degree <- igraph::degree(corGraph)
  print("Node degree:")
  print(degree)

  ## number of components

  numberOfComponents <- igraph::components(corGraph)
  numberOfComponents <- length(numberOfComponents$csize[numberOfComponents$csize > 1])

  ## all metrics

  networkMetrics <- data.frame(correlation_type = correlation_type,
                               EdgeDensity = EdgeDensity,
                               diameter = diameter,
                               averagePathLength = averagePathLength,
                               numberOfComponents = numberOfComponents,
                               averageClusteringCoefficient = averageClusteringCoefficient)

  rownames(node_label) <- node_label[, 1]
  node_label <- node_label[vars, 2]

  ### Variance results (node pie chart) ####

  if (!is.null(variance_results)) {
    rownames(variance_results) <- variance_results[, 1]
    if (!is.null(variance_results) && all(vars %in% rownames(variance_results))) {
      ps.vars <- variance_results[vars, variance_type]
      ps.vars[which(ps.vars < 0)] <- 0
      ps.vars[which(ps.vars > 1)] <- 1
    } else {
      stop("Not all variables of interest are present in the phylogenetic signal dataframe")
    }
  } else {
    ps.vars <- rep(NULL, length(correlation.matrix[, 1]))
  }

  ### Plot ####

  if(degree_as_node_size){
    nodeSize <- 5 + degree * 2
    nodeSize <- ifelse(nodeSize > 10, 10, nodeSize)
    nodeSize <- ifelse(nodeSize < 5, 5, nodeSize)
  } else {
    nodeSize <- 8
  }

  p <- qgraph::qgraph(correlation.matrix,
                      layout = layout,
                      vsize = nodeSize,
                      vsize2 = nodeSize,
                      esize = 10 * max(abs(correlation.matrix)),
                      palette = "pastel",
                      negDashed = T,
                      borders = T,
                      legend = F,
                      vTrans = 180,
                      fade = F,
                      aspect = T,
                      legend.cex = 0.25,
                      edge.labels = edge_label,
                      # label.cex = 2 * length(vars)/length(vars),
                      labels = node_label,
                      label.cex = label_size,
                      label.scale = F,
                      node.label.offset = c(0.5, -3),
                      pie = ps.vars, pieBorder = 1,
                      DoNotPlot = F,
                      groups = group_variables[, 2],
                      threshold = correlation_threshold)

  if(plot_metrics){
    graphics::text(x = 0.5, y = -1, labels = paste0("ED = ", round(networkMetrics$EdgeDensity, 2), "\n",
                                                    "NM = ", round(networkMetrics$numberOfComponents, 2), "\n",
                                                    "AC = ", round(networkMetrics$averageClusteringCoefficient, 2), "\n",
                                                    "|r|mean = ", round(networkMetrics$averagePathLength, 2), "\n",
                                                    "|r|max = ", round(networkMetrics$diameter, 2)),
                   adj = 0, cex = network_metrics_text_size)
  }
  return(p)
}
