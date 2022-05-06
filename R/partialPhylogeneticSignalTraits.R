#' Retrieve partial phylogenetic signal for each one of the traits contained in a vector once a variable or a group of variables are considered as fixed effects
#'
#' @param VARIABLES (character) Names of the variables. They must be contained in dataset.
#' @param PREDICTORS (character) Names of the predictors. They must be contained in dataset.
#' @param PHYLOGENY (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param DATASET (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param MODEL.SPECIFICATIONS (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @param FORCERUN (logical) If false, models already run are not runned again.
#'
#' @return
#' @export
#'
#' @examples
partialPhylogeneticSignalTraits <- function (VARIABLES, PREDICTORS, PHYLOGENY, DATASET, MODEL.SPECIFICATIONS = NULL,
                                             FORCERUN = T)
{
  partialPhylogeneticSignalResults <- list()
  partialPhylogeneticSignalResults$phylogenetic.signal.results <- data.frame()
  partialPhylogeneticSignalResults$models.diagnostics <- data.frame()
  partialPhylogeneticSignalResults$individual.models.results <- list()
  uni_mdls.str <- data.frame(resp_var = VARIABLES)
  uni_mdls.str$type <- paste0("partial_uni_", uni_mdls.str$resp_var)
  uni_mdls.str$n_respVars <- 1
  uni_mdls.str$pred_var <- ""
  uni_mdls.str$pred_vars <- paste0(PREDICTORS, collapse = ", ")
  uni_mdls.str$fix.frml <- paste0(uni_mdls.str$resp_var, " ~ 1")
  for (predictor in PREDICTORS) {
    uni_mdls.str$fix.frml <- paste0(uni_mdls.str$fix.frml,
                                    " + ", predictor)
  }
  uni_mdls.str$ran.frml <- "~ animal"
  uni_mdls.str$NP_ran.frml <- ""
  if (file.exists(paste0(outputs.dir, "/models_outputs/partialPhylogeneticSignalResults",
                         PREDICTORS, ".RData")) && isFALSE(FORCERUN)) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/partialPhylogeneticSignalResults",
                       PREDICTORS, ".RData"))
  }
  for (model in uni_mdls.str$type) {
    if (!model %in% names(partialPhylogeneticSignalResults$individual.models.results) |
        FORCERUN) {
      print(paste0("Running partial phylo. signal model: ",
                   model))
      model.descr <- uni_mdls.str %>% dplyr::filter(type ==
                                                      model)
      mdl.rslts <- computePartialPhylogeneticSignal(variable = model.descr$resp_var,
                                                    predictors = PREDICTORS, dataset = DATASET, phylogeny = PHYLOGENY,
                                                    model.specifications = MODEL.SPECIFICATIONS)
      partialPhylogeneticSignalResults$phylogenetic.signal.results <- rbind(partialPhylogeneticSignalResults$phylogenetic.signal.results,
                                                                            mdl.rslts$phyloSignal)
      partialPhylogeneticSignalResults$models.diagnostics <- rbind(partialPhylogeneticSignalResults$models.diagnostics,
                                                                   mdl.rslts$model.diagnostics)
      partialPhylogeneticSignalResults$individual.models.results[[model]] <- mdl.rslts
    }
  }
  print("Model structure used:")
  print(uni_mdls.str)
  print("Phylogenetic signal results:")
  print(partialPhylogeneticSignalResults$phylogenetic.signal.results)
  save(partialPhylogeneticSignalResults, file = paste0(outputs.dir,
                                                       "/models_outputs/partialPhylogeneticSignalResults", PREDICTORS,
                                                       ".RData"))
  print(paste0(outputs.dir, "/models_outputs/partialPhylogeneticSignalResults",
               PREDICTORS, ".RData"))
  return(partialPhylogeneticSignalResults)
}
