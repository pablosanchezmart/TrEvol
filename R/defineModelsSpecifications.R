#' Define mcmcglmm model specifications
#'
#' @param n_itterations (integer) Number of iterations of  Markov Chain Montecarlo.
#' @param burning (integer) Number of burning iterations.
#' @param thinning (integer) Number of thinning iterations.
#' @param uniresponse_prior (list) Parameters for the random effects of uniresponse models.
#' @param multiresponse_prior (list) Parameters for the random effects of multiresponse models
#'
#' @return List of mcmcglmm model specifications.
#' @export
#'
#' @examples
defineModelsSpecifications <- function(n_itterations = 100, burning = 10, thinning = 2,
                                          # One response variable, phylogeny and residuals
                                          uniresponse_prior = list(
                                            R = list(V = 1, nu = 0.002),
                                            G = list(G1 = list(V = 1, nu = 0.002))
                                          ),
                                          # Two response variables, phylogeny and residuals
                                          multiresponse_prior = list(
                                            R=list(V=diag(2)/2,nu=2),
                                            G=list(G1=list(V=diag(2)/2,nu=2))
                                          )
) {

  models.specifications.list <- list("number_interations" = n_itterations,
                                     "burning_iterations" = burning,
                                     "thinning_iterations" = thinning,
                                     "uniresponse_prior" = uniresponse_prior,
                                     "multiresponse_prior" = multiresponse_prior)
  return(models.specifications.list)
}
