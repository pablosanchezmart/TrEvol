computePartialCorrelations <- function(variable1, variable2, predictors = NULL, dataset, phylogeny, model.specifications = NULL) {

  if(is.null(predictors)){
    stop("Predictor needed. Define it as a character or a vector in the predictor argument.")
  }


  modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(variable1, variable2, predictors))

  # formula
  fix.frml <- paste0("cbind(", variable1, ", ", variable2, ") ~ trait-1")
  for(predictor in predictors){
    fix.frml <- paste0(fix.frml, " + trait:", predictor)
  }

  if (is.null(model.specifications)) {
    print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
    model.specifications <- defineModelsSpecifications()
  }
  # model
  mdl <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                          random = ~ us(trait):animal, rcov = ~us(trait):units,
                          data= modellingData$dta, pedigree = modellingData$phylo,
                          family = c("gaussian", "gaussian"),
                          prior = model.specifications$multiresponse_prior,
                          nitt = model.specifications$number_interations,
                          burnin = model.specifications$burning_iterations,
                          thin = model.specifications$thinning_iterations,
                          verbose = F)
  mdl$name <- fix.frml

  model.diagnostics <- diagnoseModels(model = mdl)

  # correlations calculation

  CVphylo <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".animal")]
  CVres <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".units")]

  Vphylo1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".animal")]
  Vres1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".units")]

  Vphylo2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".animal")]
  Vres2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".units")]

  # fixed effects variance
  n <- length(mdl$VCV[, 1])
  vmVarF1 <- numeric(n)
  vmVarF2 <- numeric(n)
  # separate fixed effects variance explained for each response variable
  sol_1 <- mdl$Sol[, stringr::str_detect(colnames(mdl$X), variable1)]
  X_1 <- mdl$X[, stringr::str_detect(colnames(mdl$X), variable1)]
  for(i in 1:n){
    Var <- stats::var(as.vector(sol_1[i,] %*% t(as.matrix(X_1))))
    vmVarF1[i] <- Var
  }

  sol_2 <- mdl$Sol[, stringr::str_detect(colnames(mdl$X), variable2)]
  X_2 <- mdl$X[, stringr::str_detect(colnames(mdl$X), variable2)]
  for(i in 1:n){
    Var <- stats::var(as.vector(sol_2[i,] %*% t(as.matrix(X_2))))
    vmVarF2[i] <- Var
  }


  # Total correlation (pearson correlation)

  totalCor.t1.t2  <- (CVphylo + CVres) /
    sqrt( ((Vphylo1 + Vres1 + vmVarF1) * (Vphylo2 + Vres2 + vmVarF2)) )

  correlationsChains <- data.frame("totalCorrelation" = totalCor.t1.t2)

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(totalCor.t1.t2)
  pd <- as.numeric(pd$pd)
  totalCor.t1.t2_pval <- 2*(1 - pd)

  # Partial relative phylogenetic correlation (amount of phylogenetic covariation considering total variance per trait)

  relativePhyloCor.t1.t2 <- CVphylo /
    sqrt( ((Vphylo1 + Vres1 + vmVarF1) * (Vphylo2 + Vres2 + vmVarF2)) )

  correlationsChains$parcialPhyloCorrelation <- relativePhyloCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(relativePhyloCor.t1.t2)
  pd <- as.numeric(pd$pd)
  relativePhyloCor.t1.t2_pval <- 2*(1 - pd)

  # Partial phylogenetic correlation (amount of phylogenetic covariation related to phylogenetic variation)

  phyloCor.t1.t2  <- CVphylo /
    sqrt( ((Vphylo1) * (Vphylo2)) )

  correlationsChains$PhyloCorrelation <- phyloCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(phyloCor.t1.t2)
  pd <- as.numeric(pd$pd)
  phyloCor.t1.t2_pval <- 2*(1 - pd)

  # Partial relative convergence correlation (amount of non phylogenetic covariation considering total variance per trait)

  relativeResCor.t1.t2 <- CVres /
    sqrt( (Vphylo1 + Vres1 + vmVarF1) * (Vphylo2 + Vres2 + vmVarF2) )

  correlationsChains$ParcialResCorrelation <- relativeResCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(relativeResCor.t1.t2)
  pd <- as.numeric(pd$pd)
  relativeResCor.t1.t2_pval <- 2*(1 - pd)

  # Convergence correlation (amount of residual non-phylogenetic covariation)

  resCor.t1.t2  <- CVres /
    sqrt( ((Vres1) * (Vres2)) )

  correlationsChains$ResCorrelation <- resCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(resCor.t1.t2)
  pd <- as.numeric(pd$pd)
  resCor.t1.t2_pval <- 2*(1 - pd)

  ## Results

  partialCorrelationsResults <- list()

  partialCorrelationsResults$corrrelationsSummary <- data.frame("Variable1" = variable1,
                                                                "Variable2" = variable2,
                                                                "Model" = fix.frml,
                                                                "N" = length(modellingData$dta$animal),
                                                                "Predictors" = paste0(predictors, collapse = ", "),

                                                                "Total_cor" = mean(totalCor.t1.t2),
                                                                "Phylogenetic_cor" = mean(phyloCor.t1.t2),
                                                                "Relative_phylogenetic_cor" = mean(relativePhyloCor.t1.t2),
                                                                "Convergent_cor" = mean(resCor.t1.t2),
                                                                "Relative_convergent_cor" = mean(relativeResCor.t1.t2),

                                                                "SD_Partial_total_cor" = stats::sd(totalCor.t1.t2),
                                                                "SD_phylogenetic_cor" = stats::sd(phyloCor.t1.t2),
                                                                "SD_Relative_phylogenetic_cor" = stats::sd(relativePhyloCor.t1.t2),
                                                                "SD_Convergent_cor" = stats::sd(resCor.t1.t2),
                                                                "SD_Relative_convergent_cor" = stats::sd(relativeResCor.t1.t2),

                                                                "Pvalue_Total_cor" = totalCor.t1.t2_pval,
                                                                "Pvalue_Phylogenetic_cor" = phyloCor.t1.t2_pval,
                                                                "Pvalue_Relative_phylogenetic_cor" = relativePhyloCor.t1.t2_pval,
                                                                "Pvalue_Convergent_cor" = resCor.t1.t2_pval,
                                                                "Pvalue_Relative_convergent_cor" = relativeResCor.t1.t2_pval
  )

  # Correlations
  partialCorrelationsResults$correlationsChains <- correlationsChains

  # Model
  partialCorrelationsResults$model <- mdl
  # Model diagnostics
  partialCorrelationsResults$model.diagnostics <- model.diagnostics
  return(partialCorrelationsResults)
}
