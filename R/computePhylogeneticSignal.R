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
#' @importFrom rlang .data
computePhylogeneticSignal <- function(variable, dataset, phylogeny, model.specifications = NULL,
                                       forceRun = T)
{
  modellingData <- completePhyloData(phylogeny = phylogeny,
                                     dataset = dataset, traits = variable)
  fix.frml <- paste0(variable, " ~ 1")
  if (!file.exists(paste0(.data$outputs.dir, "/models_outputs/phylo_signal.RData")) |
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
    load(file = paste0(.data$outputs.dir, "/models_outputs/phylo_signal.RData"),
         envir = .GlobalEnv)
    if(length(phylo.signal.results$phyloSignal[,1])){
      mdl <- phylo.signal.results$model
    } else{
      mdl <- phylo.signal.results$individual.models.results[[which(!is.na(stringr::str_extract(names(phylo.signal.results$individual.models.results),
                                                                                               variable)))]]$model
    }

  }
  model.diagnostics <- diagnoseModels(model = mdl)
  var.df <- modellingData$dta[, variable]
  names(var.df) <- modellingData$dta$animal
  k <- as.numeric(phytools::phylosig(modellingData$phylo, var.df,
                                     method = "K"))
  lambda <- as.numeric(phytools::phylosig(modellingData$phylo,
                                          var.df, method = "lambda")[1])
  phylo.signal.results <- list()
  phylo.signal.results$model <- mdl
  phylo.signal.results$wlambda.distr <- mdl$VCV[, "animal"]/(mdl$VCV[,
                                                                     "animal"] + mdl$VCV[, "units"])
  wlambda <- mean(phylo.signal.results$wlambda.distr)
  phylo.signal.results$resVar.distr <- mdl$VCV[, "units"]/(mdl$VCV[,
                                                                   "animal"] + mdl$VCV[, "units"])
  residualVariance <- mean(phylo.signal.results$resVar.distr)
  phylo.signal.results$phyloSignal <- data.frame(Variable = variable, N = length(modellingData$dta$animal),
                                                 Model = fix.frml, K = k, Lambda = lambda, Wlambda = wlambda,
                                                 Non_Phylogenetic_variance = residualVariance)
  phylo.signal.results$model.diagnostics <- model.diagnostics

  save(phylo.signal.results, file = paste0(.data$outputs.dir, "/models_outputs/phylo_signal.RData"))
  print(paste0(.data$outputs.dir, "/models_outputs/phylo_signal.RData"))

  return(phylo.signal.results)
}
