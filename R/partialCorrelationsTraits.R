#' Retrieve correlations for each one of the combinations among traits contained in a vector
#'
#' @param VARIABLES  (character) Names of the variables. They must be contained in dataset.
#' @param PREDICTORS (character) Names of the predictors. They must be contained in dataset.
#' @param PHYLOGENY (phylo) Phylogeny with tip labels contained in dataset$animal.
#' @param DATASET (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param MODEL.SPECIFICATIONS (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @param FORCERUN (logical) If false, models already run are not runned again.
#'
#' @return
#' @export
#'
#' @examples
partialCorrelationsTraits <- function (VARIABLES, PREDICTORS, PHYLOGENY, DATASET, MODEL.SPECIFICATIONS = NULL,
                                      FORCERUN = T) {

  # results object
  partialCorrelationsResults <- list()
  partialCorrelationsResults$correlation.results <- data.frame()
  partialCorrelationsResults$models.diagnostics <- data.frame()
  partialCorrelationsResults$individual.models.results <- list()


  # models structure

  multi_mdls.str <- expand.grid(VARIABLES, VARIABLES) # all possible pairwise combinations between variable
  names(multi_mdls.str) <- c("Var1", "Var2")

  multi_mdls.str$type <- paste0("bi_", multi_mdls.str$Var1, "_", multi_mdls.str$Var2)
  multi_mdls.str$n_respVars <- 2

  # This is needed in order to avoid models with the same response variables but in different order (which are equivalent models)
  multi_mdls.str$resp_var <- NA
  for(i in 1:length(multi_mdls.str$Var1)){
    variables <- c(as.character(multi_mdls.str[i, "Var1"]), as.character(multi_mdls.str[i, "Var2"]))
    variables <- sort(variables)
    var1 <- as.character(variables[[1]])
    var2 <- as.character(variables[[2]])
    multi_mdls.str$resp_var[i] <- paste0(var1, ", ", var2)
    multi_mdls.str$resp_var1[i] <- var1
    multi_mdls.str$resp_var2[i] <- var2
  }

  multi_mdls.str$pred_var <- ""
  multi_mdls.str$fix.frml <- paste0("cbind(", multi_mdls.str$resp_var, ") ~ trait-1")

  for(predictor in PREDICTORS){
    multi_mdls.str$fix.frml <- paste0(multi_mdls.str$fix.frml, " + trait:", predictor)
  }
  multi_mdls.str$ran.frml <- "~ us(trait):animal"
  multi_mdls.str$NP_ran.frml <- ""

  multi_mdls.str <- multi_mdls.str %>%
    dplyr::filter(!Var1 == Var2) %>%
    dplyr::filter(!duplicated(resp_var)) %>%
    dplyr::select(type, n_respVars, resp_var, resp_var1, resp_var2, pred_var, fix.frml, ran.frml, NP_ran.frml)


  # lad previous results, if exist
  if (file.exists(paste0(outputs.dir, "/models_outputs/partialCorrelationsResults", PREDICTORS, ".RData")) && isFALSE(FORCERUN)) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/partialCorrelationsResults", PREDICTORS, ".RData"))
  }

  # run models and extract results
  for (model in multi_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(partialCorrelationsResults$individual.models.results) | FORCERUN) {
      print(paste0("Running partial correlations model: ", model))
      model.descr <- multi_mdls.str %>%
      dplyr::filter(type == model)

      mdl.rslts <- computePartialCorrelations(variable1 = model.descr$resp_var1, variable2 = model.descr$resp_var2, predictors = PREDICTORS, dataset = DATASET, phylogeny = PHYLOGENY,
                                              model.specifications = MODEL.SPECIFICATIONS)
      partialCorrelationsResults$correlation.results <- rbind(partialCorrelationsResults$correlation.results,
                                                       mdl.rslts$corrrelationsSummary)
      partialCorrelationsResults$models.diagnostics <- rbind(partialCorrelationsResults$models.diagnostics,
                                                      mdl.rslts$model.diagnostics)
      partialCorrelationsResults$individual.models.results[[model]] <- mdl.rslts
    }
  }

  print("Model structure used:")
  print(multi_mdls.str)
  print("Phylogenetic signal results:")
  print(partialCorrelationsResults$phylogenetic.signal.results)

  # save results

  assign(paste0("partialCorrelationsResults", PREDICTORS), partialCorrelationsResults)

  save(mget(paste0("partialCorrelationsResults", PREDICTORS)), file = paste0(outputs.dir, "/models_outputs/partialCorrelationsResults", PREDICTORS, ".RData"))
  print(paste0(outputs.dir, "/models_outputs/partialCorrelationsResults", PREDICTORS, ".RData"))

  return(partialCorrelationsResults)
}
