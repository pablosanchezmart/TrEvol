phylogeneticSignalTraits <- function (VARIABLES, PREDICTORS, PHYLOGENY, DATASET, MODEL.SPECIFICATIONS = NULL,
                                      FORCERUN = F) {

  # results object
  partialPhylogeneticSignalResults <- list()
  partialPhylogeneticSignalResults$phylogenetic.signal.results <- data.frame()
  partialPhylogeneticSignalResults$models.diagnostics <- data.frame()
  partialPhylogeneticSignalResults$individual.models.results <- list()

  # prepare models structure
  uni_mdls.str <- data.frame(resp_var = VARIABLES)
  uni_mdls.str$type <- paste0("partial_uni_", uni_mdls.str$resp_var)
  uni_mdls.str$n_respVars <- 1
  uni_mdls.str$pred_var <- ""
  uni_mdls.str$pred_vars <- paste0(PREDICTORS, collapse = ", ")
  for(predictor in PREDICTORS){
    res_uni_mdls.str$fix.frml <- paste0(res_uni_mdls.str$fix.frml, " + ", predictor)
  }
  uni_mdls.str$ran.frml <- "~ animal"
  uni_mdls.str$NP_ran.frml <- ""

  # lad previous results, if exist
  if (file.exists(paste0(outputs.dir, "/models_outputs/partialPhylogeneticSignalResults", PREDICTORS, ".RData"))) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/partialPhylogeneticSignalResults", PREDICTORS, ".RData"))
  }

  # run models and extract results
  for (model in uni_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(partialPhylogeneticSignalResults$individual.models.results) | FORCERUN) {
      print(paste0("Running phylo. signal model: ", model))
      model.descr <- uni_mdls.str %>%
        dplyr::filter(type == model)
      mdl.rslts <- computePhylogeneticSignal(variable = model.descr$resp_var, predictors = PREDICTORS, dataset = DATASET, phylogeny = PHYLOGENY, model.specifications = MODEL.SPECIFICATIONS)

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

  # save results
  save(partialPhylogeneticSignalResults, file = paste0(outputs.dir, "/models_outputs/partialPhylogeneticSignalResults", PREDICTORS, ".RData"))
  print(paste0(outputs.dir, "/models_outputs/partialPhylogeneticSignalResults", PREDICTORS, ".RData"))

  return(partialPhylogeneticSignalResults)
}
