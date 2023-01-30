#' Plot variance-covariance matrix (phylogenetic signal and relative phylogentic and convergent correlations)
#'
#' @param correlations (data frame) correlations output of correlationsTraits$correlations.rslts
#' @param phylogenetic.signal (data frame) phylogenetic signal as reported by phylogeneticSignalTraits$phylogenetic.signal.rslts
#' @param order_vars (character) order of the variables to plot.
#' @param labels (data frame) name of the variables and their corresponding labels to be ploted.
#' @param triangular (logical) if true, all the information is plotted together in an efficient way in a triangular matrix.
#' @param triangularFillType (pie or "circle) Ow to represent proportional correlation"
#'
#' @return
#' @export
#'
#' @examples
plotVcv <- function (correlations, phylogenetic.signal, order_vars = NULL, labels = NULL, triangular = T, triangularFillType = "pie") {

  # select significant correlations
  correlations[which(correlations[, "Pvalue_Total_coordination"] > 0.05 &
                       correlations[, "Pvalue_Total_coordinated_phylogenetic_conservatism"] >
                       0.05), c("Pvalue_Total_coordination", "Total_coordinated_phylogenetic_conservatism",
                                "Total_coordinated_radiation")] <- 0

  # Calculate proportional correlations (relative to total correlation)

  correlations$cor_proportion <- abs(correlations[, "Total_coordinated_phylogenetic_conservatism"])/(abs(correlations[, "Total_coordinated_radiation"]) + (abs(correlations[, "Total_coordinated_phylogenetic_conservatism"])))
  correlations$cor_proportion[which(is.na(correlations$cor_proportion))] <- 0

  # Prepare matrices
  correlations <- correlations %>% dplyr::select(Trait_1,
                                                 Trait_2, cor_proportion, Total_coordination)
  vars <- unique(c(correlations$Trait_1, correlations$Trait_2))
  proportion.correlation.matrix <- matrix(ncol = length(vars),
                                          nrow = length(vars), 0)
  colnames(proportion.correlation.matrix) <- vars
  rownames(proportion.correlation.matrix) <- vars
  total.correlation.matrix <- proportion.correlation.matrix
  for (i in 1:length(correlations$Trait_1)) {
    var1 <- correlations$Trait_1[i]
    var2 <- correlations$Trait_2[i]
    proportion.correlation.matrix[var1, var2] <- correlations[i,
                                                              "cor_proportion"]
    proportion.correlation.matrix[var2, var1] <- correlations[i,
                                                              "cor_proportion"]
    proportion.correlation.matrix[var1, var1] <- phylogenetic.signal[phylogenetic.signal$Trait ==
                                                                       var1, "Total_phylogenetic_conservatism"]
    proportion.correlation.matrix[var2, var2] <- phylogenetic.signal[phylogenetic.signal$Trait ==
                                                                       var2, "Total_phylogenetic_conservatism"]
    total.correlation.matrix[var1, var2] <- correlations[i,
                                                         "Total_coordination"]
    total.correlation.matrix[var2, var1] <- correlations[i,
                                                         "Total_coordination"]
    total.correlation.matrix[var1, var1] <- 1
    total.correlation.matrix[var2, var2] <- 1
  }

  # order traits

    if (!is.null(order_vars)) {
    proportion.correlation.matrix <- proportion.correlation.matrix[order_vars,
                                                                   order_vars]
    total.correlation.matrix <- total.correlation.matrix[order_vars,
                                                         order_vars]
    }

  # Change labels

  if (!is.null(labels)) {
    rownames(labels) <- labels[, 1]
    labels <- labels[rownames(total.correlation.matrix), ] %>%
      dplyr::pull(name)
    rownames(proportion.correlation.matrix) <- labels
    colnames(proportion.correlation.matrix) <- labels
    rownames(total.correlation.matrix) <- labels
    colnames(total.correlation.matrix) <- labels
  }

  # Plot

  if (isTRUE(triangular)) {
    p <- TrEvol:::plotVcvTriangular(corr = total.correlation.matrix,
                                    corrProp = proportion.correlation.matrix, fillType = triangularFillType)
  } else {
    mixed.matrix <- proportion.correlation.matrix
    mixed.matrix[upper.tri(mixed.matrix, diag = F)] <- total.correlation.matrix[upper.tri(total.correlation.matrix, diag = F)]
    mixed.matrix[lower.tri(mixed.matrix, diag = T)] <- proportion.correlation.matrix[lower.tri(proportion.correlation.matrix, diag = T)]
    p <- corrplot::corrplot.mixed(mixed.matrix, lower = "pie",
                                  upper = "circle", diag = "l", tl.pos = "lt", tl.col = "black")
  }

  return(p)
}

