#' Plot variance-covariance matrix (phylogenetic signal and relative phylogentic and convergent correlations)
#'
#' @param correlations (data frame) correlations output of correlationsTraits$correlations.rslts
#' @param phylogenetic.signal (data frame) phylogenetic signal as reported by phylogeneticSignalTraits$phylogenetic.signal.rslts
#' @param order_vars (character) order of the variables to plot.
#'
#' @return
#' @export
#'
#' @examples
plotVcv <- function (correlations, phylogenetic.signal, order_vars = NULL) {
  correlations[which(correlations[, "Pvalue_Total_cor"] > 0.05 & correlations[, "Pvalue_Relative_phylogenetic_cor"] > 0.05),
               c("Total_cor", "Relative_phylogenetic_cor")] <- 0
  correlations <- correlations[, c("Variable1", "Variable2", "Total_cor", "Relative_phylogenetic_cor")]

  correlations$cor_proportion <- abs(correlations[, "Total_cor"])/(abs(correlations[, "Total_cor"]) + (abs(correlations[, "Relative_phylogenetic_cor"])))
  correlations$cor_proportion[which(is.na(correlations$cor_proportion))] <- 0
  correlations <- correlations %>% dplyr::select(Variable1, Variable2, cor_proportion)
  vars <- unique(c(correlations$Variable1, correlations$Variable2))
  correlation.matrix <- matrix(ncol = length(vars), nrow = length(vars), 0)
  colnames(correlation.matrix) <- vars
  rownames(correlation.matrix) <- vars

  for (i in 1:length(correlations$Variable1)) {
    var1 <- correlations$Variable1[i]
    var2 <- correlations$Variable2[i]
    correlation.matrix[var1, var2] <- correlations[i, "cor_proportion"]
    correlation.matrix[var2, var1] <- correlations[i, "cor_proportion"]
    correlation.matrix[var1, var1] <- phylogenetic.signal[phylogenetic.signal$Variable == var1, "Wlambda"]
    correlation.matrix[var2, var2] <- phylogenetic.signal[phylogenetic.signal$Variable == var2, "Wlambda"]
  }
  if (!is.null(order_vars)) {
    correlation.matrix <- correlation.matrix[order_vars,
                                             order_vars]
  }
  p <- corrplot::corrplot(correlation.matrix, method = "pie", type = "lower", diag = T, cl.pos = "n", tl.col = "black")
  return(p)
}
