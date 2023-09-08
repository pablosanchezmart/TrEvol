#' Simulate a list of photogenically correlated and non-photogenically correlated traits given a phylogeny
#'
#'Returns a list with a data frame containing the simulated traits and the variance-covariance matrix used. The number of traits simulated doubles the
#'number of rows in the supplied variance-covariance matrix, which needs to be symmetric (i.e., same number of rows and columns).
#'This is because two sets of traits are simulated under the same variance covariance. In the first set of traits,
#'their variances and covariances are phylogenetically conserved, while the second set of traits variances
#'and covariances are phylogenetically independent.
#'
#' @param phylogeny (*phylo*). Phylogeny under which simulate the phylogenetically correlated traits (i.e., phylogenetically conserved traits). If not specified,
#' the function generates a phylogeny using the function simulate_bm_model from the  castor R package. When not specified, the number of observations
#' need to be specified in the number_observations argument.
#' @param number_observations (*integer*). Number of observations to simulate if phylogeny is not provided. When the phylogeny is provided, the number of terminal taxa of
#' the phylogeny is used, and then, this argument is ignored.
#' @param trait.names (*character*). Names of the traits to simulate following the variance covariance order (i.e., first name correspond to the first row and column).
#' @param vcv_matrix (*matrix*). Variance covariance matrix setting the covariance among traits. The number of rows and colums needs to match the lenght as trait.names.
#'
#' @return A list containing a data frame with phylogenetically and non-phylogenetically correlated traits, the variance-covariance matrix and the phylogeny used.
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate phylogenetically correlated and phylogenetically independent traits using default parameters
#' simulatedTraits.data <- simulateDataSet()
#' }
simulateDataSet <- function(phylogeny = NULL,
                            number_observations = 100,
                            trait_names = c("G1_trait1", "G1_trait2", "G1_env", "G2_trait1", "G2_trait2", "G2_env"),
                            vcv_matrix = matrix(c(1, 0.9, 0.8, 0, 0.1, 0.2,
                                                 0.9, 1, 0.8, 0, 0.1, 0.2,
                                                 0.8, 0.8, 1, 0, 0.1, 0.2,
                                                 0, 0, 0, 1, 0.9, 0.8,
                                                 0.1, 0.1, 0.1, 0.9, 1, 0.8,
                                                 0.2, 0.2, 0.2, 0.8, 0.8, 1), ncol = 6))  {

  # simulate phylogeny if it is not provided
  if(is.null(phylogeny)){
    phylogeny <- phytools::pbtree(n = number_observations)
  }

  # Check whether variance covariance is provided and if it's symmetric

  if(!is.matrix(vcv_matrix)){
    stop("Please provide a symmetric variance covariance matrix object in the vcv_matrix. is.matrix(vcv_matrix) needs to be TRUE.")
  } else{
    if(isSymmetric.matrix(vcv_matrix)){
      stop("Please provide a symmetric variance covariance matrix object in the vcv_matrix. isSymmetric.matrix(vcv_matrix) needs to be TRUE.")
    }
  }

  # Check whether the trait.names exist and if it is the same length as the number of rows and columns in the variance-covariance matrix

  if(length(trait.names) != dim(vcv_matrix)){
    stop("Please provide a vector with trait names. The length of the vector needs to be equal to the number of columns and rows of the variance covariance matrix.")
  }

  # Create phylogenetically correlated and independent names
  BM_trait.names <- paste0("phylo_", trait.names)
  nonBM_trait.names <- paste0("nonPhylo_", trait.names)

  # Name the matrix
  diffMat <- vcv_matrix
  colnames(diffMat) <- BM_trait.names
  rownames(diffMat) <- BM_trait.names

  # Simulate phylogenetically correlated traits (under a Brownian motion model of evolution)
  BM.df = as.data.frame(castor::simulate_bm_model(phylogeny, diffusivity=diffMat)$tip_states)

  # Name simulated traits
  colnames(BM.df) <- BM_trait.names

  # Add the terminal taxa names as a column and row names to match phylogeny and data
  BM.df <- cbind("animal" = phylogeny$tip.label, BM.df)
  rownames(BM.df) <- BM.df$animal

  # Simulate non-phylogenetically correlated traits
  nonBM.df <- faux::rnorm_multi(n = length(phylogeny$tip.label),
                          mu = c(0, 0, 0, 0, 0, 0),
                          sd = c(1, 1, 1, 1, 1, 1),
                          r = vcv_matrix,
                          varnames = nonBM_trait.names,
                          empirical = FALSE)

  # Add terminal taxa names
  nonBM.df <- cbind("animal" = phylogeny$tip.label, nonBM.df)
  rownames(nonBM.df) <- nonBM.df$animal

  # Merge phylogenetically and non-phylogenetically correlated traits
  df <- merge(BM.df, nonBM.df, by = "animal")

  # Print the variance-covariance matrix used
  print(diffMat)

  # Report results
  rslts <- list()

  rslts$data <- df
  rslts$vcv_matrix <- diffMat
  rslts$phylogeny <- phylogeny
  return(rslts)
}
