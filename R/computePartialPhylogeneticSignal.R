computePartialPhylogeneticSignal <- function(variable, predictors = NULL, dataset, phylogeny, model.specifications = NULL) {

  if(is.null(predictors)){
    stop("Predictor needed. Define it as a character or a vector in the predictor argument.")
  }


  modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(variable, predictors))

  # formula
  fix.frml <- paste0(variable, "~ 1 ")
  for(predictor in predictors){
    fix.frml <- paste0(fix.frml, " + ", predictor)
  }

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

  # Residual phylogenetic signal
  n <- length(mdl$VCV[, 1])
  vmVarF <- numeric(n)

  # fixed effect explained variance (from Kakagawa 2013)
  for(i in 1:n){
    Var <- var(as.vector(mdl$Sol[i,] %*% t(as.matrix(mdl$X))))
    vmVarF[i] <- Var
  }

  # partial phylogenetic signal calculation
  wlambda.distr <- (mdl$VCV[, "animal"])/(vmVarF + mdl$VCV[, "animal"] + mdl$VCV[, "units"])

  resNonPhylogeneticVar <- (mdl$VCV[, "units"])/(vmVarF + mdl$VCV[, "animal"] + mdl$VCV[, "units"])

  # results
  partialPhylogeneticSignalResults <- list()
  partialPhylogeneticSignalResults$phyloSignal <- data.frame("Variable" = response,
                                                      "N" = length(modellingData$dta$animal),
                                                      "Model" = fix.frml,
                                                      "Partial_Wlambda" = mean(wlambda.distr),
                                                      "Partial_non_phylogenetic_variance" = mean(resNonPhylogeneticVar))
  partialPhylogeneticSignalResults$wlambda.distr <- wlambda.distr
  partialPhylogeneticSignalResults$resNonPhylogeneticVar <- resNonPhylogeneticVar
  partialPhylogeneticSignalResults$model <- mdl
  partialPhylogeneticSignalResults$model.diagnostics <- model.diagnostics

  return(partialPhylogeneticSignalResults)
}
