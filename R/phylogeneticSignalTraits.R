#' Retrieve phylogenetic signal for each one of the traits contained in a vector
#'
#' @param VARIABLES (character) Names of the variables. They must be contained in dataset.
#' @param PHYLOGENY (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param DATASET (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param MODEL.SPECIFICATIONS (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @param FORCERUN (logical) If false, models already run are not runned again.
#'
#' @return
#' @export
#'
#' @examples
phylogeneticSignalTraits <- function (VARIABLES, PHYLOGENY, DATASET, MODEL.SPECIFICATIONS = NULL,
                                      FORCERUN = F) {

  # results object
  phylogeneticSignalResults <- list()
  phylogeneticSignalResults$phylogenetic.signal.results <- data.frame()
  phylogeneticSignalResults$models.diagnostics <- data.frame()
  phylogeneticSignalResults$individual.models.results <- list()

  # prepare models structure
  uni_mdls.str <- data.frame(resp_var = VARIABLES)
  uni_mdls.str$type <- paste0("uni_", uni_mdls.str$resp_var)
  uni_mdls.str$n_respVars <- 1
  uni_mdls.str$pred_var <- ""
  uni_mdls.str$fix.frml <- paste0(uni_mdls.str$resp_var, " ~ 1")
  uni_mdls.str$ran.frml <- "~ animal"

  # lad previous results, if exist
  if (file.exists(paste0(outputs.dir, "/models_outputs/phylogeneticSignalResults.RData")) && isFALSE(FORCERUN)) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/phylogeneticSignalResults.RData"))
  }

  # run models and extract results
  for (model in uni_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(phylogeneticSignalResults$individual.models.results) | FORCERUN) {
      print(paste0("Running phylo. signal model: ", model))
      model.descr <- uni_mdls.str %>%
        dplyr::filter(type == model)
      mdl.rslts <- computePhylogeneticSignal(variable = model.descr$resp_var, dataset = DATASET, phylogeny = PHYLOGENY, model.specifications = MODEL.SPECIFICATIONS)

      phylogeneticSignalResults$phylogenetic.signal.results <- rbind(phylogeneticSignalResults$phylogenetic.signal.results,
                                                                     mdl.rslts$phyloSignal)
      phylogeneticSignalResults$models.diagnostics <- rbind(phylogeneticSignalResults$models.diagnostics,
                                                            mdl.rslts$model.diagnostics)
      phylogeneticSignalResults$individual.models.results[[model]] <- mdl.rslts
    }
  }

  print("Model structure used:")
  print(uni_mdls.str)
  print("Phylogenetic signal results:")
  print(phylogeneticSignalResults$phylogenetic.signal.results)

  # save results
  save(phylogeneticSignalResults, file = paste0(outputs.dir, "/models_outputs/phylogeneticSignalResults.RData"))
  print(paste0(outputs.dir, "/models_outputs/phylogeneticSignalResults.RData"))

  return(phylogeneticSignalResults)
}
