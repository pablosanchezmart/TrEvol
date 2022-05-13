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
  correlations[which(correlations[, "Pvalue_Total_cor"] > 0.05 & correlations[, "Pvalue_Relative_phylogenetic_cor"] > 0.05), c("Total_cor", "Relative_phylogenetic_cor", "Relative_convergent_cor")] <- 0

  correlations$cor_proportion <- abs(correlations[, "Relative_phylogenetic_cor"]) /
    (abs(correlations[, "Relative_convergent_cor"]) + (abs(correlations[, "Relative_phylogenetic_cor"])))

  correlations$cor_proportion[which(is.na(correlations$cor_proportion))] <- 0
  correlations <- correlations %>% dplyr::select(Variable1,
                                                 Variable2, cor_proportion, Total_cor)
  vars <- unique(c(correlations$Variable1, correlations$Variable2))

  proportion.correlation.matrix <- matrix(ncol = length(vars), nrow = length(vars),
                                          0)
  colnames(proportion.correlation.matrix) <- vars
  rownames(proportion.correlation.matrix) <- vars
  total.correlation.matrix <- proportion.correlation.matrix

  for (i in 1:length(correlations$Variable1)) {
    var1 <- correlations$Variable1[i]
    var2 <- correlations$Variable2[i]
    proportion.correlation.matrix[var1, var2] <- correlations[i, "cor_proportion"]
    proportion.correlation.matrix[var2, var1] <- correlations[i, "cor_proportion"]
    proportion.correlation.matrix[var1, var1] <- phylogenetic.signal[phylogenetic.signal$Variable ==
                                                                       var1, "Wlambda"]
    proportion.correlation.matrix[var2, var2] <- phylogenetic.signal[phylogenetic.signal$Variable ==
                                                                       var2, "Wlambda"]

    total.correlation.matrix[var1, var2] <- correlations[i, "Total_cor"]
    total.correlation.matrix[var2, var1] <- correlations[i, "Total_cor"]
    total.correlation.matrix[var1, var1] <- 1
    total.correlation.matrix[var2, var2] <- 1
  }

  if (!is.null(order_vars)) {
    proportion.correlation.matrix <- proportion.correlation.matrix[order_vars,
                                                                   order_vars]

    total.correlation.matrix <- total.correlation.matrix[order_vars,
                                                         order_vars]
  }

  if(!is.null(labels)){
    rownames(labels) <- labels[, 1]
    labels <- labels[rownames(total.correlation.matrix), ] %>% dplyr::pull(name)

    rownames(proportion.correlation.matrix) <- labels
    colnames(proportion.correlation.matrix) <- labels

    rownames(total.correlation.matrix) <- labels
    colnames(total.correlation.matrix) <- labels
  }

  if(isTRUE(triangular)){
    p <- plotVcvTriangular(corr = total.correlation.matrix, corrProp = proportion.correlation.matrix, fillType = triangularFillType)
  } else{
    mixed.matrix <- proportion.correlation.matrix
    mixed.matrix[upper.tri(mixed.matrix, diag = F)] <- total.correlation.matrix[upper.tri(total.correlation.matrix, diag = F)]
    mixed.matrix[lower.tri(mixed.matrix, diag = T)] <- proportion.correlation.matrix[lower.tri(proportion.correlation.matrix, diag = T)]

    p <- corrplot::corrplot.mixed(mixed.matrix, lower = 'pie', upper = 'circle', diag = "l", tl.pos = "lt", tl.col = "black")
  }


  return(p)
}

