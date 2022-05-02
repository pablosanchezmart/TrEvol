#' Compute phylogenetic signal using univariate phylogenetic mixed models.
#'
#' @param variable (character) Name of the variable. It must be contained in datset.
#' @param dataset (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param model.specifications Mcmcglmm models specifications. See defineModelsSpecification.
#' @param forceRun (logical) If false, models already run are not runned again.
#'
#' @return List containing univariate model results, diagnostics and phylogenetic signal results summary. Phlogenetic signal using Blomberg's K and Pagel's lambda are also reported (using phytools phylosig() fun)
#' @export
#'
#' @examples
#' @importFrom rlang .data
computePhylogeneticSignal <- function(variable, dataset, phylogeny, model.specifications = NULL, forceRun = T){


  # Complete data
  modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = variable)
  # Model formula
  fix.frml <- paste0(variable, " ~ 1")

  # Run models if previous results do not exist
  if(!file.exists(paste0(outputs.dir, "/models_outputs/phyloSignalResults.RData")) | forceRun){

    # Model specifications
    if(is.null(model.specifications)){
      print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
      model.specifications <- defineModelsSpecifications()
      # print(model.specifications)
    }

    # Model
    mdl <- MCMCglmm::MCMCglmm(stats::as.formula(fix.frml), random = ~ animal,
                              family = "gaussian",
                              prior = model.specifications$uniresponse_prior,
                              data = modellingData$dta,
                              pedigree = modellingData$phylo,
                              nitt = model.specifications$number_interations,
                              burnin = model.specifications$burning_iterations,
                              thin = model.specifications$thinning_iterations,
                              verbose=F)
    mdl$name <- fix.frml

  } else {
    print("loanding previous results")
    load(file = paste0(.data$outputs.dir, "/models_outputs/phyloSignalResults.RData"),  envir = .GlobalEnv)
    mdl <- .data$phylo.signal.results$individual.models.results[[which(!is.na(stringr::str_extract(names(phylo.signal.results$individual.models.results), variable)))]]$model
  }

  # Bayesian model diagnostics check
  model.diagnostics <- diagnoseModels(model = mdl)

  # Phylo. signal using phytools

  var.df <- modellingData$dta[, variable]
  names(var.df) <- modellingData$dta$animal

  k <- as.numeric(phytools::phylosig(modellingData$phylo, var.df, method = "K"))
  lambda <- as.numeric(phytools::phylosig(modellingData$phylo, var.df, method = "lambda")[1])

  ## Results gathering
  rslts <- list()
  rslts$model <- mdl

  # Weighted lambda from mcmcglmm results
  rslts$wlambda.distr <- mdl$VCV[,"animal"]/ (mdl$VCV[, "animal"] + mdl$VCV[, "units"])
  wlambda <- mean(rslts$wlambda.distr)

  # Weighted residual variance from mcmcglmm results
  rslts$resVar.distr <- mdl$VCV[, "units"]/ (mdl$VCV[, "animal"] + mdl$VCV[, "units"])
  residualVariance <- mean(rslts$resVar.distr)

  # Phylogenetic signal data frame
  rslts$phyloSignal <- data.frame("Variable" = variable,
                                  "N" = length(modellingData$dta$animal),
                                  "Model" = fix.frml,
                                  "K" = k,
                                  "Lambda" = lambda,
                                  "Wlambda" = wlambda,
                                  "Non_Phylogenetic_variance" = residualVariance)

  # Model diagnostics
  rslts$model.diagnostics <- model.diagnostics

  return(rslts)
}
