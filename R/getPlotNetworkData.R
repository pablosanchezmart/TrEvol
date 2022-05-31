#' Get data to plot networks
#'
#' @param covariance.results (data frame) correlations output of computeCovariancePartition()$covarianceResults
#' @param variance.results (data frame) phylogenetic signal as reported by computeVariancePartition()$varianceResults
#' @param variance.type (character) variable name of the phylogenetic signal as reported in varianceResults
#' @param covariance.type (character) variable name of the correlation as reported in correlations
#' @param only.significant  (logical) if TRUE, only significant correlations are plotted
#'
#' @return
#' @export
#'
#' @examples
getPlotNetworkData <- function(covariance.results, variance.results = NULL, variance.type = NULL, covariance.type = "Total_coordination",
                               only.significant = T){


  data.list <- list()
  # correlation matrix

  # set non significant values to zero if only.significant = T
  if (only.significant) {
    covariance.results[which(covariance.results[, paste0("Pvalue_", covariance.type)] > 0.05), covariance.type] <- 0
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

  data.list$Edges <- correlation.matrix

  # Variance partition (node size)
  if (!is.null(variance.results)) {
    rownames(variance.results) <- variance.results[, 1]
    if (!is.null(variance.results) && all(vars %in% rownames(variance.results))) {
      ps.vars <- variance.results[vars, variance.type]
      ps.vars[which(ps.vars < 0)] <- 0
      ps.vars[which(ps.vars > 1)] <- 1
      names(ps.vars) <- vars
    } else {
      stop("Not all variables of interest are present in the phylogenetic signal dataframe")
    }
  } else {
    ps.vars <- rep(1, length(correlation.matrix[, 1]))
  }

  data.list$NodeSize <- ps.vars

  # Node coordinates (x, y and z)

  pcoanalisis <- ecodist::pco(correlation.matrix[lower.tri(correlation.matrix)])
  data.list$coordinates <- pcoanalisis$vectors[, 1:3]
  rownames(data.list$coordinates) <- vars
  colnames(data.list$coordinates) <- c("x", "y", "z")

  return(data.list)
}
