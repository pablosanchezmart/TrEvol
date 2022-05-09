plotCorrelationMatrices <- function(correlations,
                                    phylogenetic.signal,
                            correlation.type1 = "Proportional_phylo_cor_mean",
                            correlation.type2 = "Proportional_resid_cor_mean",
                            varsOrder = NULL){

  # Plot only significant correlations
  correlations[which(correlations[ , paste0("Total_cor", "_pvalue")] > 0.05 & correlations[ , paste0("Total_cor", "_pvalue")] > 0.05), c("Total_cor", "Phylogenetic_cor")] <- 0

  correlations <- correlations[, c("Variable1", "Variable2", "Total_cor", "Phylogenetic_cor")]

  # Phylogenetic correlation proportion over total correlation
  correlations$cor_proportion <- abs(correlations[, "Total_cor"]) / (abs(correlations[, "Total_cor"]) + (abs(correlations[, "Phylogenetic_cor"])))
  correlations$cor_proportion[which(is.na(correlations$cor_proportion))] <- 0
  correlations <- correlations %>% dplyr::select(Variable1, Variable2, cor_proportion)

  # Correlation matrix

  vars <- unique(c(correlations$Variable1, correlations$Variable2))
  correlation.matrix <- matrix(ncol = length(vars), nrow = length(vars), 0)
  colnames(correlation.matrix) <- vars
  rownames(correlation.matrix) <- vars

  for(i in 1:length(correlations$Variable1)){
    var1 <- correlations$Variable1[i]
    var2 <- correlations$Variable2[i]
    correlation.matrix[var1, var2] <- correlations[i, "cor_proportion"]
    correlation.matrix[var2, var1] <- correlations[i, "cor_proportion"]
    correlation.matrix[var1, var1] <- phylogenetic.signal[phylogenetic.signal$Variable == var1, "Wlambda"]
    correlation.matrix[var2, var2] <- phylogenetic.signal[phylogenetic.signal$Variable == var2, "Wlambda"]
  }

  if(!is.null(varsOrder)){
    correlation.matrix <- correlation.matrix[varsOrder, varsOrder]
  }

  p <- corrplot::corrplot(correlation.matrix, lower = "circle", upper="pie", diag = T, cl.pos = "n", tl.col = "black")
  return(p)
}
