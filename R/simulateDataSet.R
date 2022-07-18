#' Simulate a list of phylogenetically correlated and non-phylogenetically correlated traits given a phylogeny
#'
#' @param phylogeny phylogeny under wich simulate phylogenetically correlated traits
#' @param nObs Number of observations to simulate if phylogeny is not provided. If provided, the number of terminal taxa of the phylogeny is used.
#' @param trait.names Names of the traits to simulate
#' @param vcvMatrix Variance covariance matrix setting the covariance among traits. Ncol have to be the same lenght as trait.names.
#'
#' @return A data frame with phylogenetically and non-phylogenetically correlated traits.
#' @export
#'
#' @examples
simulateDataSet <- function(phylogeny = NULL, nObs = 100, trait.names = c("G1_trait1", "G1_trait2", "G1_envPred", "G2_trait1", "G2_trait2", "G2_envPred"),
                                                               vcvMatrix = matrix(c(1, 0.9, 0.8, 0, 0.1, 0.2,
                                                                                    0.9, 1, 0.8, 0, 0.1, 0.2,
                                                                                    0.8, 0.8, 1, 0, 0.1, 0.2,
                                                                                    0, 0, 0, 1, 0.9, 0.8,
                                                                                    0.1, 0.1, 0.1, 0.9, 1, 0.8,
                                                                                    0.2, 0.2, 0.2, 0.8, 0.8, 1), ncol = 6))  {
  # simulate phylogeny if it is not provided
  if(is.null(phylogeny)){
    phylogeny <- phytools::pbtree(n = nObs)
  }

  BM_trait.names <- paste0("phylo_", trait.names)
  nonBM_trait.names <- paste0("nonPhylo_", trait.names)

  diffMat <- vcvMatrix
  colnames(diffMat) <- BM_trait.names
  rownames(diffMat) <- BM_trait.names

  BM.df = as.data.frame(castor::simulate_bm_model(phylogeny, diffusivity=diffMat)$tip_states)
  colnames(BM.df) <- BM_trait.names
  BM.df <- cbind("animal" = phylogeny$tip.label, BM.df)

  rownames(BM.df) <- BM.df$animal

  nonBM.df <- faux::rnorm_multi(n = length(phylogeny$tip.label),
                          mu = c(0, 0, 0, 0, 0, 0),
                          sd = c(1, 1, 1, 1, 1, 1),
                          r = vcvMatrix,
                          varnames = nonBM_trait.names,
                          empirical = FALSE)


  nonBM.df <- cbind("animal" = phylogeny$tip.label, nonBM.df)

  rownames(nonBM.df) <- nonBM.df$animal

  df <- merge(BM.df, nonBM.df, by = "animal")
  return(df)
}
