#' Compute correlations using bivariate phylogenetic mixed models for a given pair of variables.
#'
#' @param variable1 (character) Name of the variable. It must be contained in datset.
#' @param variable2 (character) Name of the variable. It must be contained in datset.
#' @param dataset (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param model.specifications (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @param forceRun (logical) If false, models already run are not runned again.
#'
#' @return
#' @export
#'
#' @examples
computeCorrelations <- function(variable1, variable2, dataset, phylogeny,  model.specifications = NULL, forceRun = F) {

  # Model specifications
  if(is.null(model.specifications)){
    print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
    model.specifications <- defineModelsSpecifications()
    # print(model.specifications)
  }
  # Complete data
  modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(variable1, variable2))

  # Model
  fix.frml <- paste0("cbind(", variable1, ", ", variable2, ") ~ trait-1")

  # Run models if previous results do not exist
  if(!file.exists(paste0(outputs.dir, "/models_outputs/correlationResults.RData")) | forceRun){

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

  } else {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/correlationsResults.RData"))

    if(length(correlationsResults$corrrelationsSummary[, 1]) > 1){
    mdl <- correlationsResults$individual.models.results[[which(!is.na(stringr::str_extract(names(correlationsResults$individual.models.results), variable1)) & !is.na(stringr::str_extract(names(correlationsResults$individual.models.results), variable2)))]]$model
    } else{
      mdl <- correlationsResults$model
    }

}
  # Bayesian model diagnostics check
  model.diagnostics <- diagnoseModels(mdl)

  ## Correlations
  CVphylo <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".animal")]
  CVres <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".units")]

  Vphylo1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".animal")]
  Vres1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".units")]

  Vphylo2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".animal")]
  Vres2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".units")]

  # Total correlation (pearson correlation)
  totalCor.t1.t2  <- (CVphylo + CVres) /
    sqrt( (Vphylo1 + Vres1) * (Vphylo2 + Vres2) )

  correlationsChains <- data.frame("totalCorrelation" = totalCor.t1.t2)

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(totalCor.t1.t2)
  pd <- as.numeric(pd$pd)
  totalCor.t1.t2_pval <- 2*(1 - pd)

  # Partial phylogenetic correlation (amount of phylogenetic covariation considering total variance per trait)
  relativePhyloCor.t1.t2 <- CVphylo /
    sqrt( (Vphylo1 + Vres1) * (Vphylo2 + Vres2) )

  correlationsChains$parcialPhyloCorrelation <- relativePhyloCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(relativePhyloCor.t1.t2)
  pd <- as.numeric(pd$pd)
  relativePhyloCor.t1.t2_pval <- 2*(1 - pd)

  # Phylogenetic correlation (amount of phylogenetic covariation related to phylogenetic variation)
  phyloCor.t1.t2  <- CVphylo /
    sqrt( (Vphylo1 * Vphylo2) )

  correlationsChains$PhyloCorrelation <- phyloCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(phyloCor.t1.t2)
  pd <- as.numeric(pd$pd)
  phyloCor.t1.t2_pval <- 2*(1 - pd)

  # Partial residual correlation (amount of residual covariation considering total variance per trait)
  relativeResCor.t1.t2 <- CVres /
    sqrt( (Vphylo1 + Vres1) * (Vphylo2 + Vres2) )

  correlationsChains$ParcialResCorrelation <- relativeResCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(relativeResCor.t1.t2)
  pd <- as.numeric(pd$pd)
  relativeResCor.t1.t2_pval <- 2*(1 - pd)

  # Residual correlation (amount of residual covariation related to residual trait variation)

  resCor.t1.t2  <- CVres /
    sqrt( (Vres1 * Vres2) )

  correlationsChains$ResCorrelation <- resCor.t1.t2

  # pvalue related to correlation being different than 0
  pd <- bayestestR::p_direction(resCor.t1.t2)
  pd <- as.numeric(pd$pd)
  resCor.t1.t2_pval <- 2*(1 - pd)

  ## Results gathering
  correlationsResults <- list()

  t1.t2.cor.results <- data.frame("Variable1" = variable1,
                                  "Variable2" = variable2,
                                  "Model" = fix.frml,
                                  "N" = length(modellingData$dta$animal),

                                  "Total_cor" = mean(totalCor.t1.t2),
                                  "Phylogenetic_cor" = mean(phyloCor.t1.t2),
                                  "Proportional_phylogenetic_cor" = mean(relativePhyloCor.t1.t2),
                                  "Convergent_cor" = mean(resCor.t1.t2),
                                  "Proportional_convergent_cor" = mean(relativeResCor.t1.t2),

                                  "SD_Total_cor" = stats::sd(totalCor.t1.t2),
                                  "SD_Phylogenetic_cor" = stats::sd(phyloCor.t1.t2),
                                  "SD_Proportional_phylogenetic_cor" = stats::sd(relativePhyloCor.t1.t2),
                                  "SD_Convergent_cor_sd" = stats::sd(resCor.t1.t2),
                                  "SD_Proportional_convergent_cor" = stats::sd(relativeResCor.t1.t2),

                                  "Pvalue_Total_cor" = totalCor.t1.t2_pval,
                                  "Pvalue_Phylogenetic_cor" = phyloCor.t1.t2_pval,
                                  "Pvalue_Proportional_phylogenetic_cor" = relativePhyloCor.t1.t2_pval,
                                  "Pvalue_Convergent_cor" = resCor.t1.t2_pval,
                                  "Pvalue_Proportional_convergent_cor" = relativeResCor.t1.t2_pval
  )

  ## Results gathering

  # Model
  correlationsResults$model <- mdl
  # Correlations
  correlationsResults$correlations.distr <- correlationsChains
  correlationsResults$corrrelationsSummary <- t1.t2.cor.results
  # Model diagnostics
  correlationsResults$model.diagnostics <- model.diagnostics
  return(correlationsResults)
}
