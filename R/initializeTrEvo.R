#' Create and define results and outputs folders structure.
#'
#' @param folderName (character) Folder name.
#'
#' @return
#' @export
#'
#' @examples
initializeTrEvo <- function(folderName = "simulations"){

  print("loading workspace ...")
  options(scipen = 999, digits = 3)
  wd <- getwd()
  outputs.dir <- paste0(wd, "/outputs/outputs_", folderName)
  results.dir <- paste0(wd, "/results/results_", folderName)
  dir.create("outputs")
  dir.create(outputs.dir)
  dir.create(paste0(outputs.dir, "/models_outputs"))
  dir.create(paste0(outputs.dir, "/phylogenetic_variance_covariance"))
  dir.create(paste0(outputs.dir, "/predictions"))
  dir.create("results")
  dir.create(results.dir)
  dir.create(paste0(results.dir, "/figures"))
  dir.create(paste0(results.dir, "/figures/phylogeneticSignal_plots/"))
  dir.create(paste0(results.dir, "/figures/VCV_networks/"))
  dir.create(paste0(results.dir, "/figures/conditional_VCV_networks/"))
  dir.create(paste0(results.dir, "/figures/example/"))
  dir.create(paste0(results.dir, "/figures/predictions/"))
  dir.create(paste0(results.dir, "/tables"))
  envList <- list("outputs.dir", "results.dir")

  assign("outputs.dir", outputs.dir, envir = .GlobalEnv)
  assign("results.dir", results.dir, envir = .GlobalEnv)
}
