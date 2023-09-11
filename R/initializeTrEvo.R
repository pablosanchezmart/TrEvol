#' Create and define results and outputs folders structure.
#'
#'Create folders within the working directory where results and outputs of TrEvol functions will be stored if "save" arguments are set to TRUE.
#'Note that this function creates directories in your working directory and will assign their location in the objects "output.dir" and "results.dir".
#' @param folder_name (*character*) Name that will be given to the folder within which the results and outputs subfolders will be created.
#'
#' @return
#' @examples
#' \dontrun{
#' # Create outputs and results folders.
#'
#' # Simulate example data
#' initializeTrEvo(folder_name = "simulations")
#'
#' # Outputs directory
#' outputs.dir
#'
#' Results directory
#' results.dir
#'
#' }
#'
#' @examples
initializeTrEvo <- function(folder_name = "simulations"){

  print("loading workspace ...")
  options(scipen = 999, digits = 3)

  # Get user working directory
  wd <- getwd()

  # Create outputs and results dir if not yet created, otherwise, ignored

  outputs.dir <- paste0(wd, "/outputs/outputs_", folder_name)
  results.dir <- paste0(wd, "/results/results_", folder_name)

  dir.create("outputs", showWarnings = F)
  dir.create(outputs.dir, showWarnings = F)
  dir.create(paste0(outputs.dir, "/models_outputs"), showWarnings = F)
  dir.create(paste0(outputs.dir, "/phylogenetic_variance_covariance"), showWarnings = F)
  dir.create(paste0(outputs.dir, "/predictions"), showWarnings = F)

  dir.create("results", showWarnings = F)
  dir.create(results.dir, showWarnings = F)
  dir.create(paste0(results.dir, "/figures"), showWarnings = F)
  dir.create(paste0(results.dir, "/figures/VCV_networks/"), showWarnings = F)
  dir.create(paste0(results.dir, "/figures/scatterplots/"), showWarnings = F)
  dir.create(paste0(results.dir, "/figures/predictions/"), showWarnings = F)
  dir.create(paste0(results.dir, "/tables"), showWarnings = F)

  # Assign directories to objects to further use
  envList <- list("outputs.dir", "results.dir")
  assign("outputs.dir", outputs.dir, envir = .GlobalEnv)
  assign("results.dir", results.dir, envir = .GlobalEnv)
}
