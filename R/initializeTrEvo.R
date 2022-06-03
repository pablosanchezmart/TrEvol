#' Create and define results and outputs folders structure.
#'
#' @param folder.name (character) Folder name.
#'
#' @return
#' @export
#'
#' @examples
initializeTrEvo <- function(folder.name = "simulations"){

  print("loading workspace ...")
  options(scipen = 999, digits = 3)
  wd <- getwd()
  outputs.dir <- paste0(wd, "/outputs/outputs_", folder.name)
  results.dir <- paste0(wd, "/results/results_", folder.name)
  dir.create("outputs", showWarnings = F)
  dir.create(outputs.dir, showWarnings = F)
  dir.create(paste0(outputs.dir, "/models_outputs"), showWarnings = F)
  dir.create(paste0(outputs.dir, "/predictions"), showWarnings = F)
  dir.create("results", showWarnings = F)
  dir.create(results.dir, showWarnings = F)
  dir.create(paste0(results.dir, "/figures"), showWarnings = F)
  dir.create(paste0(results.dir, "/figures/VCV_networks/"), showWarnings = F)
  dir.create(paste0(results.dir, "/figures/scatterplots/"), showWarnings = F)
  dir.create(paste0(results.dir, "/figures/predictions/"), showWarnings = F)
  dir.create(paste0(results.dir, "/tables"), showWarnings = F)

  envList <- list("outputs.dir", "results.dir")
  assign("outputs.dir", outputs.dir, envir = .GlobalEnv)
  assign("results.dir", results.dir, envir = .GlobalEnv)
}
