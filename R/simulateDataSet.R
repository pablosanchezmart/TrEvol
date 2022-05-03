#' Simulate a list of phylogenetically correlated and non-phylogenetically correlated traits given a phylogeny
#'
#' @param phylogeny phylogeny under wich simulate phylogenetically correlated traits
#' @param nObs Number of observations to simulate if phylogeny is not provided. If provided, the number of terminal taxa of the phylogeny is used.
#' @param traitNames Names of the traits to simulate
#' @param vcvMatrix Variance covariance matrix setting the covariance among traits. Ncol have to be the same lenght as traitNames.
#'
#' @return A data frame with phylogenetically and non-phylogenetically correlated traits.
#' @export
#'
#' @examples
simulateDataSet <- function(phylogeny = NA, nObs = 100, traitNames = c("HC_1", "HC_2", "HC_predictor", "LC_1",  "LC_2",   "LC_predictor"),
                                                               vcvMatrix = matrix(c(1, 0.9, 0.8, 0, 0.1, 0.2,
                                                                                    0.9, 1, 0.8, 0, 0.1, 0.2,
                                                                                    0.8, 0.8, 1, 0, 0.1, 0.2,
                                                                                    0, 0, 0, 1, 0.9, 0.8,
                                                                                    0.1, 0.1, 0.1, 0.9, 1, 0.8,
                                                                                    0.2, 0.2, 0.2, 0.8, 0.8, 1), ncol = 6))  {
  # simulate phylogeny if it is not provided
  if(is.na(phylogeny)){
    phylogeny <- phytools::pbtree(n = nObs)
  }

  BM_traitNames <- paste0("BM_", traitNames)
  nonBM_traitNames <- paste0("nonBM_", traitNames)

  diffMat <- vcvMatrix
  colnames(diffMat) <- BM_traitNames
  rownames(diffMat) <- BM_traitNames

  BM.df = as.data.frame(castor::simulate_bm_model(phylogeny, diffusivity=diffMat)$tip_states)
  colnames(BM.df) <- BM_traitNames
  BM.df <- cbind("animal" = phylogeny$tip.label, BM.df)

  rownames(BM.df) <- BM.df$animal

  nonBM.df <- faux::rnorm_multi(n = length(phylogeny$tip.label),
                          mu = c(0, 0, 0, 0, 0, 0),
                          sd = c(1, 1, 1, 1, 1, 1),
                          r = diffMat,
                          varnames = nonBM_traitNames,
                          empirical = FALSE)


  nonBM.df <- cbind("animal" = phylogeny$tip.label, nonBM.df)

  rownames(nonBM.df) <- nonBM.df$animal

  df <- merge(BM.df, nonBM.df, by = "animal")
  return(df)
}
