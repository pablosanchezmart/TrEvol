#' Compute phylogenetic signal using univariate phylogenetic mixed models for a given variable.
#'
#' @param variable (character) Name of the variable. It must be contained in dataset.
#' @param dataset (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param model.specifications (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @param forceRun (logical) If false, models already run are not runned again.
#'
#' @return List containing univariate model results, diagnostics and phylogenetic signal results summary. Phlogenetic signal using Blomberg's K and Pagel's lambda are also reported (using phytools phylosig() fun)
#' @export
computePhylogeneticSignal <- function(variable, dataset, phylogeny, model.specifications = NULL,
                                       forceRun = T)
{
  modellingData <- completePhyloData(phylogeny = phylogeny,
                                     dataset = dataset, traits = variable)
  fix.frml <- paste0(variable, " ~ 1")
  if (!file.exists(paste0(outputs.dir, "/models_outputs/phylogeneticSignalResults.RData")) |
      forceRun) {
    if (is.null(model.specifications)) {
      print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
      model.specifications <- defineModelsSpecifications()
    }
    mdl <- MCMCglmm::MCMCglmm(stats::as.formula(fix.frml),
                              random = ~animal, family = "gaussian", prior = model.specifications$uniresponse_prior,
                              data = modellingData$dta, pedigree = modellingData$phylo,
                              nitt = model.specifications$number_interations, burnin = model.specifications$burning_iterations,
                              thin = model.specifications$thinning_iterations,
                              verbose = F)
    mdl$name <- fix.frml
  }
  else {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/phylogeneticSignalResults.RData"),
         envir = .GlobalEnv)
    if(length(phylogeneticSignalResults$phyloSignal[, 1]) > 1){
      mdl <- phylogeneticSignalResults$individual.models.results[[which(!is.na(stringr::str_extract(names(phylogeneticSignalResults$individual.models.results),
                                                                                               variable)))]]$model
    } else{
      mdl <- phylogeneticSignalResults$model
    }

  }
  model.diagnostics <- diagnoseModels(model = mdl)
  var.df <- modellingData$dta[, variable]
  names(var.df) <- modellingData$dta$animal
  k <- as.numeric(phytools::phylosig(modellingData$phylo, var.df,
                                     method = "K"))
  lambda <- as.numeric(phytools::phylosig(modellingData$phylo,
                                          var.df, method = "lambda")[1])
  phylogeneticSignalResults <- list()
  phylogeneticSignalResults$model <- mdl
  phylogeneticSignalResults$wlambda.distr <- mdl$VCV[, "animal"]/(mdl$VCV[,
                                                                     "animal"] + mdl$VCV[, "units"])
  wlambda <- mean(phylogeneticSignalResults$wlambda.distr)
  phylogeneticSignalResults$resVar.distr <- mdl$VCV[, "units"]/(mdl$VCV[,
                                                                   "animal"] + mdl$VCV[, "units"])
  residualVariance <- mean(phylogeneticSignalResults$resVar.distr)
  phylogeneticSignalResults$phyloSignal <- data.frame(Variable = variable,
                                                 N = length(modellingData$dta$animal),
                                                 Model = fix.frml,
                                                 K = k,
                                                 Lambda = lambda,
                                                 Wlambda = wlambda,
                                                 Non_Phylogenetic_variance = residualVariance)
  phylogeneticSignalResults$model.diagnostics <- model.diagnostics

  return(phylogeneticSignalResults)
}
