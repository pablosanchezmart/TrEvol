computePhylogeneticSignal <- function(variable, dataset, phylogeny, model.specifications = NULL) {

  modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = variable)
  fix.frml <- paste0(variable, " ~ 1")

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

  model.diagnostics <- diagnoseModels(model = mdl)
  var.df <- modellingData$dta[, variable]
  names(var.df) <- modellingData$dta$animal
  k <- as.numeric(phytools::phylosig(modellingData$phylo, var.df,
                                     method = "K"))
  lambda <- as.numeric(phytools::phylosig(modellingData$phylo,
                                          var.df, method = "lambda")[1])

  wlambda.distr <- mdl$VCV[, "animal"]/(mdl$VCV[, "animal"] + mdl$VCV[, "units"])
  resVar.distr <- mdl$VCV[, "units"]/(mdl$VCV[, "animal"] + mdl$VCV[, "units"])

  # reuslts
  phylogeneticSignalResults <- list()
  phylogeneticSignalResults$phyloSignal <- data.frame("Variable" = variable,
                                                      "N" = length(modellingData$dta$animal),
                                                      "Model" = fix.frml,
                                                      "K" = k,
                                                      "Lambda" = lambda,
                                                      "Wlambda" = mean(wlambda.distr),
                                                      "Non_Phylogenetic_variance" = mean(resVar.distr))

  phylogeneticSignalResults$wlambda.distr <- wlambda.distr
  phylogeneticSignalResults$resVar.distr <- resVar.distr
  phylogeneticSignalResults$model <- mdl
  phylogeneticSignalResults$model.diagnostics <- model.diagnostics

  return(phylogeneticSignalResults)
}
