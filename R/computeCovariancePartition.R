#' Compute covarince partition
#'
#' @param trait1
#' @param trait2
#' @param environmentalVariables
#' @param dataset
#' @param phylogeny
#' @param model.specifications
#'
#' @return
#' @export
#'
#' @examples

computeCovariancePartition <- function(traits = c("BM_HC_1", "BM_HC_2"), environmentalVariables = "BM_HC_predictor", dataset, phylogeny,
                                       model.specifications = modelSpecifications) {


  # results object
  traitsCovariancePartitionResults <- list()
  traitsCovariancePartitionResults$covarianceResults <- data.frame()
  traitsCovariancePartitionResults$models.diagnostics <- data.frame()
  traitsCovariancePartitionResults$individual.models.results <- list()

  # models structure

  multi_mdls.str <- expand.grid(traits, traits) # all possible pairwise combinations between variable
  names(multi_mdls.str) <- c("Var1", "Var2")

  multi_mdls.str$type <- paste0("bi_", multi_mdls.str$Var1, "_", multi_mdls.str$Var2)

  # This is needed in order to avoid models with the same response variables but in different order (which are equivalent models)
  multi_mdls.str$resp_var <- NA
  for(i in 1:length(multi_mdls.str$Var1)){
    traits <- c(as.character(multi_mdls.str[i, "Var1"]), as.character(multi_mdls.str[i, "Var2"]))
    traits <- sort(traits)
    var1 <- as.character(traits[[1]])
    var2 <- as.character(traits[[2]])
    multi_mdls.str$resp_var[i] <- paste0(var1, ", ", var2)
    multi_mdls.str$resp_var1[i] <- var1
    multi_mdls.str$resp_var2[i] <- var2
  }

  multi_mdls.str$fix.frml <- paste0("cbind(", multi_mdls.str$resp_var, ") ~ trait-1")
  multi_mdls.str$ran.frml <- "~ us(trait):animal"

  multi_mdls.str <- multi_mdls.str %>%
    dplyr::filter(!Var1 == Var2) %>%
    dplyr::filter(!duplicated(resp_var)) %>%
    dplyr::select(type, resp_var, resp_var1, resp_var2, fix.frml, ran.frml)

  # lad previous results, if exist
  if (file.exists(paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults", environmentalVariables, ".RData")) && isFALSE(FORCERUN)) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults", environmentalVariables, ".RData"))
  }

  # run models and extract results
  for (model in multi_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(traitsCovariancePartitionResults$individual.models.results) | FORCERUN) {

      print(paste0("Running covariance calculation: ", model))
      model.descr <- multi_mdls.str %>%
        dplyr::filter(type == model)

      trait1 <- model.descr$resp_var1
      trait2 <- model.descr$resp_var2

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait1, trait2, environmentalVariables))

      # formula
      fix.frml <- paste0("cbind(", trait1, ", ", trait2, ") ~ trait-1")

      if (is.null(model.specifications)) {
        print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
        model.specifications <- defineModelsSpecifications()
      }

      ### Covariance partition without accounting for environment ####

      ## Phylogenetic model

      mdlPhylo <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                     random = ~ us(trait):animal, rcov = ~us(trait):units,
                                     data= modellingData$dta, pedigree = modellingData$phylo,
                                     family = c("gaussian", "gaussian"),
                                     prior = model.specifications$multiresponse_prior,
                                     nitt = model.specifications$number_interations,
                                     burnin = model.specifications$burning_iterations,
                                     thin = model.specifications$thinning_iterations,
                                     verbose = F)

      mdlPhylo$name <- fix.frml

      model.diagnostics <- diagnoseModels(model = mdlPhylo)

      ### Covariance partition calculation

      # total covariance
      totalCV <- mdlPhylo$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")] + mdlPhylo$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]

      # total phyogenetic covariance
      totalPhyloCV <- mdlPhylo$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")]

      # phylogenetic variance
      totalPhyloVar1 <- mdlPhylo$VCV[, paste0("trait", trait1, ":trait", trait1, ".animal")]
      totalPhyloVar2 <- mdlPhylo$VCV[, paste0("trait", trait2, ":trait", trait2, ".animal")]

      # total residual variance
      totalResidualCV <- mdlPhylo$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]

      # non-phylogenetic variance
      totalNonPhyloVar1 <- mdl$VCV[, paste0("trait", trait2, ":trait", trait2, ".units")]
      totalNonPhyloVar2 <- mdl$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")]

      ## results

      # total correlation

      totalCoordination  <- totalCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ) )

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(totalCoordination)
      pd <- as.numeric(pd$pd)
      totalCoordination_pval <- 2*(1 - pd)

      # total coordinated phylogenetic conservatism

      totalCoordinatedPhylogeneticConservatism <- totalPhyloCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(totalCoordinatedPhylogeneticConservatism)
      pd <- as.numeric(pd$pd)
      totalCoordinatedPhylogeneticConservatism_pval <- 2*(1 - pd)

      # total non-phylogenetic coordination

      totalCoordinatedRadiation  <- totalResidualCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(totalCoordinatedRadiation)
      pd <- as.numeric(pd$pd)
      totalCoordinatedRadiation_pval <- 2*(1 - pd)

      covariancePartitionResults <- list()

      covariancePartitionResults$covariancePartition <-  data.frame("Trait_1" = trait1,
                                                                    "Trait_2" = trait2,
                                                                    "N" = length(modellingData$dta$animal),

                                                                    "Total_coordination" = mean(totalCoordination),
                                                                    "Total_coordinated_phylogenetic_conservatism" = mean(totalCoordinatedPhylogeneticConservatism),
                                                                    "Total_coordinated_radiation" = mean(totalCoordinatedRadiation),

                                                                    "Pvalue_Total_coordination" = totalCoordination_pval,
                                                                    "Pvalue_Total_coordinated_phylogenetic_conservatism" = totalCoordinatedPhylogeneticConservatism_pval,
                                                                    "Pvalue_Total_coordinated_radiation" = totalCoordinatedRadiation_pval

      )
      covariancePartitionResults$covariancePartitionDistributions <- list("totalCoordination" = totalCoordination,
                                                                          "totalCoordinatedPhylogeneticConservatism" = totalCoordinatedPhylogeneticConservatism,
                                                                          "Total_coordinated_radiation" = Total_coordinated_radiation)

      covariancePartitionResults$modelPhylo <- mdlPhylo
      covariancePartitionResults$model.diagnostics <- model.diagnostics


      if(is.null(environmentalVariables)){
      # add to all traits results
      traitsCovariancePartitionResults$covarianceResults <- rbind(traitsCovariancePartitionResults$covarianceResults,
                                                              covariancePartitionResults$covariancePartition)

      traitsCovariancePartitionResults$models.diagnostics <- rbind(traitsCovariancePartitionResults$models.diagnostics,
                                                             covariancePartitionResults$model.diagnostics)
      traitsCovariancePartitionResults$individual.models.results[[model]] <- covariancePartitionResults
      }


      ### Variance partition including the environment ####

      if(!is.null(environmentalVariables)){

      for(predictor in environmentalVariables){
        fix.frml <- paste0(fix.frml, " + trait:", predictor)
      }

      mdlPhyloEnv <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                        random = ~ us(trait):animal, rcov = ~us(trait):units,
                                        data= modellingData$dta, pedigree = modellingData$phylo,
                                        family = c("gaussian", "gaussian"),
                                        prior = model.specifications$multiresponse_prior,
                                        nitt = model.specifications$number_interations,
                                        burnin = model.specifications$burning_iterations,
                                        thin = model.specifications$thinning_iterations,
                                        verbose = F)

      mdlPhyloEnv$name <- fix.frml

      model.diagnostics <- diagnoseModels(model = mdlPhyloEnv)


      ### Covariance partition calculation

      # pure phylogenetic covariance
      purePhyloCV <-  mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")]

      # avoid negative variances when purePhylo covariance is a little bit higher than totalPhylo variance
      # if(mean(purePhyloCV) > mean(totalPhyloCV)){
      #   purePhyloCV <- totalPhyloCV
      # }

      # pure residual covariance
      pureResidualCV <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]


      # total coordinated phylogenetic conservatism

      pureCoordinatedPhylogeneticConservatism <- purePhyloCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(pureCoordinatedPhylogeneticConservatism)
      pd <- as.numeric(pd$pd)
      pureCoordinatedPhylogeneticConservatism_pval <- 2*(1 - pd)

      # total coordinated phylogenetic niche conservatism

      coordinatedNichePhylogeneticConservatism <- (totalPhyloCV - purePhyloCV) / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(coordinatedNichePhylogeneticConservatism)
      pd <- as.numeric(pd$pd)
      coordinatedNichePhylogeneticConservatism_pval <- 2*(1 - pd)


      # pure environmental coordination

      pureEnvironmentalCoordination <- (totalResidualCV - pureResidualCV) / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(pureEnvironmentalCoordination)
      pd <- as.numeric(pd$pd)
      pureEnvironmentalCoordination_pval <- 2*(1 - pd)

      totalEnvironmentalCoordination <- coordinatedNichePhylogeneticConservatism + pureEnvironmentalCoordination

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(totalEnvironmentalCoordination)
      pd <- as.numeric(pd$pd)
      totalEnvironmentalCoordination_pval <- 2*(1 - pd)

      residualCoordination <- 1 - (pureCoordinatedPhylogeneticConservatism + pureCoordinatedNichePhylogeneticConservatism + pureEnvironmentalCoordination)

      # pvalue related to correlation being different than 0
      pd <- bayestestR::p_direction(residualCoordination)
      pd <- as.numeric(pd$pd)
      residualCoordination_pval <- 2*(1 - pd)

      # results

      covariancePartitionResults$covariancePartition <-  cbind(covariancePartitionResults$covariancePartition,
                                                               "Environmental_variables" = paste0(environmentalVariables, collapse = ", "),
                                                               "Pure_coordinated_phylogenetic_conservatism" = mean(pureCoordinatedPhylogeneticConservatism),
                                                               "Coordinated_niche_phylogenetic_conservatism" = mean(coordinatedNichePhylogeneticConservatism),
                                                               "Total_environmental_coordination" = mean(totalEnvironmentalCoordination),
                                                               "Pure_environmental_coordination" = mean(pureEnvironmentalCoordination),
                                                               "Residual_coordination" = mean(residualCoordination),

                                                               "Pvalue_Pure_coordinated_phylogenetic_conservatism" = pureCoordinatedPhylogeneticConservatism_pval,
                                                               "Pvalue_Coordinated_niche_phylogenetic_conservatism" = coordinatedNichePhylogeneticConservatism_pval,
                                                               "Pvalue_Total_environmental_coordination" = totalEnvironmentalCoordination_pval,
                                                               "Pvalue_Pure_environmental_coordination" = pureEnvironmentalCoordination_pval,
                                                               "residualCoordination" = residualCoordination_pval
      )

      covariancePartitionResults$covariancePartitionDistributions["pureCoordinatedPhylogeneticConservatism"] <- pureCoordinatedPhylogeneticConservatism
      covariancePartitionResults$covariancePartitionDistributions["coordinatedNichePhylogeneticConservatism"] <- coordinatedNichePhylogeneticConservatism
      covariancePartitionResults$covariancePartitionDistributions["totalEnvironmentalCoordination"] <- totalEnvironmentalCoordination
      covariancePartitionResults$covariancePartitionDistributions["pureEnvironmentalCoordination"] <- pureEnvironmentalCoordination
      covariancePartitionResults$covariancePartitionDistributions["residualCoordination"] <- residualCoordination

      covariancePartitionResults$modelPhyloEnv <- mdlPhyloEnv
      covariancePartitionResults$model.diagnostics <- model.diagnostics

      # add to all traits results
      traitsCovariancePartitionResults$covarianceResults <- rbind(traitsCovariancePartitionResults$covarianceResults,
                                                              covariancePartitionResults$covariancePartition)

      traitsCovariancePartitionResults$models.diagnostics <- rbind(traitsCovariancePartitionResults$models.diagnostics,
                                                             covariancePartitionResults$model.diagnostics)
      traitsCovariancePartitionResults$individual.models.results[[model]] <- covariancePartitionResults

      } # end evaluation if model already exists in results

    } # end bucle for all traits
  }

  print("Model structure used:")
  print(multi_mdls.str)
  print("Phylogenetic signal results:")
  print(traitsCovariancePartitionResults$covarianceResults)

  # save results

  save(list = "traitsCovariancePartitionResults", file = paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults.RData"))
  print(paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults.RData"))

  return(traitsCovariancePartitionResults)
}
