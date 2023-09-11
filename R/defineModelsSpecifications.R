#' Define *mcmcglmm* (Bayesian Phylogenetic Mixed Models) specifications
#'
#' Return an object containing the specifications that can be used as in the argument "model_specifications" of the functions
#' computeVariancePartition, computeCovariancePartition and computeVarianceCovariancePartition.
#' See "MCMCglmm" R package documentation for further information.
#'
#' @param number_iterations (*integer*) Number of iterations of  the MCMCglmm model.
#' @param burning (*integer*) Number of burning iterations. This initial iterations will be excluded from when estimating the posterior distribution. It allows achieving convergence.
#' @param thinning (*integer*) Number of thinning iterations. Number of iterations that will be discarded after a sampled one. It allows to deal with autocorrelation.
#' @param uniresponse_prior (*list*) Parameters for the random effects for uni-response models. Default is set to an inverse-Gamma distribution, which is canonical and considered non-informative.
#' @param biresponse_prior (*list*) Parameters for the random effects of bi-response models. Normally, these are the models including two traits. Default is set to an inverse-Gamma distribution, which is canonical and considered non-informative.
#' @param triresponse_prior (*list*) Parameters for the random effects of tri-response models. Normally, these are the models including two traits and one environmental variable. Default is set to an inverse-Gamma distribution, which is canonical and considered non-informative.
#'
#' @return List of mcmcglmm model specifications.
#' @export
#'
#' @examples
#' \dontrun{
#' # Define mcmcglmm model specifications
#' model_specifications <- defineModelsSpecifications(
#' number_iterations = 100,
#' burning = 10,
#' thinning = 2,
#'
#'# One response variable, phylogeny and residuals
#'uniresponse_prior = list(
#'  R = list(V = 1, nu = 0.002),
#'  G = list(G1 = list(V = 1, nu = 0.002))
#'),
#'# Two response variables, phylogeny and residuals
#'biresponse_prior = list(
#'  R=list(V=diag(2)/2,nu=2),
#'  G=list(G1=list(V=diag(2)/2,nu=2))
#'),
#'# Two response variables, phylogeny and residuals
#'triresponse_prior = list(
#'  R=list(V=diag(3)/2,nu=2),
#'  G=list(G1=list(V=diag(3)/2,nu=2)))
# )
# }
defineModelsSpecifications <- function(number_iterations = 100, burning = 10, thinning = 2,
                                          # One response variable, phylogeny and residuals
                                          uniresponse_prior = list(
                                            R = list(V = 1, nu = 0.002),
                                            G = list(G1 = list(V = 1, nu = 0.002))
                                          ),
                                       # Two response variables, phylogeny and residuals
                                       biresponse_prior = list(
                                         R=list(V=diag(2)/2,nu=2),
                                         G=list(G1=list(V=diag(2)/2,nu=2))
                                       ),
                                       # Three response variables, phylogeny and residuals
                                       triresponse_prior = list(
                                         R=list(V=diag(3)/2,nu=2),
                                         G=list(G1=list(V=diag(3)/2,nu=2))
                                       )
) {

  # Model specification list
  models.specifications.list <- list("number_interations" = number_iterations,
                                     "burning_iterations" = burning,
                                     "thinning_iterations" = thinning,
                                     "uniresponse_prior" = uniresponse_prior,
                                     "biresponse_prior" = biresponse_prior,
                                     "triresponse_prior" = triresponse_prior)
  return(models.specifications.list)
}
