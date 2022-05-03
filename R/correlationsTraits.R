#' Retrieve correlations for each one of the combinations among traits contained in a vector
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
correlationsTraits <- function(VARIABLES, PHYLOGENY, DATASET, MODEL.SPECIFICATIONS = NULL, FORCERUN = F){

  allCorrelation.rslts <- data.frame()
  # models structure

  multi_mdls.str <- expand.grid(variables.list, variables.list) # all possible pairwise combinations between variables

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
  multi_mdls.str$ran.frml <- "~ us(trait):animal"
  multi_mdls.str$NP_ran.frml <- ""
  multi_mdls.str$pred_var <- ""
  multi_mdls.str$fix.frml <- paste0("cbind(", multi_mdls.str$resp_var, ") ~ trait-1")
  multi_mdls.str$ran.frml <- "~ us(trait):animal"
  multi_mdls.str$NP_ran.frml <- ""

  multi_mdls.str <- multi_mdls.str %>%
    dplyr::filter(!Var1 == Var2) %>%
    dplyr::filter(!duplicated(resp_var)) %>%
    dplyr::select(type, n_respVars, resp_var, resp_var1, resp_var2, pred_var, fix.frml, ran.frml, NP_ran.frml)

  correlationsResults <- list()
  correlationsResults$allCorrelation.rslts <- data.frame()
  correlationsResults$allModelDiagnostics <- data.frame()
  correlationsResults$individual.rslts <- list()

  # load previous results
  if (file.exists(paste0(outputs.dir, "/models_outputs/correlationsResults.RData"))) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/correlationsResults.RData"))
  }

  # run models not included in results
  for(model in multi_mdls.str$type){
    if (!model %in% names(correlationsResults$individual.models.results) | FORCERUN) {
    print(paste0("Running correlations model: ", model))

    model.descr <- multi_mdls.str %>%
      dplyr::filter(type == model)

    mdl.rslts <- computeCorrelations(variable1 = model.descr$resp_var1, variable2 = model.descr$resp_var2, datasert = DATASET, phylogeny = PHYLOGENY,
                                     model.specifications = MODEL.SPECIFICATIONS)
    correlationsResults$individual.models.results[[model]] <-  mdl.rslts
    correlationsResults$correlation.results <- rbind(correlationsResults$correlation.results, mdl.rslts$corrrelationsSummary)
    correlationsResults$models.diagnostics <- rbind(correlationsResults$models.diagnostics, mdl.rslts$model.diagnostics)
    }
  }
  correlationsResults("Model structure used:")
  print(multi_mdls.str)
  print("Phylogenetic signal results:")
  print(allCorrelation.rslts)
  return(correlationsResults)
}
