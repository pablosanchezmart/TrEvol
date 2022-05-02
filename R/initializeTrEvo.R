#' Initialize package and set outputs and results directories.
#'
#' @param folderName (Charaacter) name of the folder where results and outputs will be stored.
#'
#' @return
#' @export
#'
#' @examples
initializeTrEvo <- function(folderName = "simulations"){

  print("loading workspace ...")

  ### GENERAL SPECIFICATIONS ----------------------------------------------------- ####

  options(scipen = 999, digits = 3)

  # set the different paths for results using a subset of the data depending on the grouping variable
  folderName <- paste0(getwd(), "/", folderName)
  print(paste0("Outputs and results will be stored in: ", folderName))

  outputs.dir <- paste0("outputs/outputs_", folderName)
  results.dir <- paste0("results/results_", folderName)

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
}


