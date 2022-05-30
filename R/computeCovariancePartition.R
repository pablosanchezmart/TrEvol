#' Covariance partition including (or not) environment
#'
#' @param traits (character) Name of the trait or list of traits. It  must be contained in the dataset.
#' @param environmental.variables (character) Names of the environmental variables They must be contained in dataset.
#' @param dataset (data frame) Dataset containing the trait of interest and a column named "animal" describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param model.specifications (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @param force.run (logical) If false, models already run are not runned again.
#' @param save (logical) If false, resulta re not saved.
#'
#' @return
#' @export
#'
#' @examples

computeCovariancePartition <- function(traits, environmental.variables = NULL, dataset, phylogeny, model.specifications = NULL, force.run = T, save = T) {


  # results object
  traitsCovariancePartitionResults <- list()
  traitsCovariancePartitionResults$covarianceResults <- data.frame()
  traitsCovariancePartitionResults$models.diagnostics <- data.frame()
  traitsCovariancePartitionResults$individual.models.results <- list()

  # models structure

  multi_mdls.str <- expand.grid(traits, traits) # all possible pairwise combinations between variable
  names(multi_mdls.str) <- c("trait1", "trait2")

  multi_mdls.str$type <- paste0("bi_", multi_mdls.str$trait1, "_", multi_mdls.str$trait2)

  # This is needed in order to avoid models with the same response variables but in different order (which are equivalent models)
  multi_mdls.str$traits <- NA
  for(i in 1:length(multi_mdls.str$trait1)){
    traits <- c(as.character(multi_mdls.str[i, "trait1"]), as.character(multi_mdls.str[i, "trait2"]))
    traits <- sort(traits)
    multi_mdls.str$traits[i] <- paste0(traits[1], ", ", traits[2])
    multi_mdls.str$trait1[i] <- as.character(traits[[1]])
    multi_mdls.str$trait2[i] <-  as.character(traits[[2]])
  }

  multi_mdls.str$fix.frml <- paste0("cbind(", multi_mdls.str$traits, ") ~ trait-1")
  multi_mdls.str$ran.frml <- "~ us(trait):animal"

  multi_mdls.str <- multi_mdls.str %>%
    dplyr::filter(!trait1 == trait2) %>%
    dplyr::filter(!duplicated(traits)) %>%
    dplyr::select(type, traits, trait1, trait2, fix.frml, ran.frml)

  # lad previous results, if exist
  if (file.exists(paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults", environmental.variables, ".RData")) && isFALSE(force.run)) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults", environmental.variables, ".RData"))
  }

  # run models and extract results
  for (model in multi_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(traitsCovariancePartitionResults$individual.models.results) | force.run) {

      print(paste0("Running covariance calculation: ", model))
      model.descr <- multi_mdls.str %>%
        dplyr::filter(type == model)

      trait1 <- as.character(model.descr$trait1)
      trait2 <- as.character(model.descr$trait2)

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait1, trait2, environmental.variables))

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
      totalNonPhyloVar1 <- mdlPhylo$VCV[, paste0("trait", trait2, ":trait", trait2, ".units")]
      totalNonPhyloVar2 <- mdlPhylo$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")]

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
                                                                          "Total_coordinated_radiation" = totalCoordinatedRadiation)

      covariancePartitionResults$modelPhylo <- mdlPhylo
      covariancePartitionResults$model.diagnostics <- model.diagnostics


      if(is.null(environmental.variables)){
        # add to all traits results
        traitsCovariancePartitionResults$covarianceResults <- rbind(traitsCovariancePartitionResults$covarianceResults,
                                                                    covariancePartitionResults$covariancePartition)

        traitsCovariancePartitionResults$models.diagnostics <- rbind(traitsCovariancePartitionResults$models.diagnostics,
                                                                     covariancePartitionResults$model.diagnostics)
        traitsCovariancePartitionResults$individual.models.results[[model]] <- covariancePartitionResults
      }


      ### Variance partition including the environment ####

      if(!is.null(environmental.variables)){

        for(predictor in environmental.variables){
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

        residualCoordination <- 1 - (pureCoordinatedPhylogeneticConservatism + coordinatedNichePhylogeneticConservatism + pureEnvironmentalCoordination)

        # pvalue related to correlation being different than 0
        pd <- bayestestR::p_direction(residualCoordination)
        pd <- as.numeric(pd$pd)
        residualCoordination_pval <- 2*(1 - pd)

        # results

        covariancePartitionResults$covariancePartition <-  cbind(covariancePartitionResults$covariancePartition,
                                                                 "Environmental_variables" = paste0(environmental.variables, collapse = ", "),
                                                                 "Pure_coordinated_phylogenetic_conservatism" = mean(pureCoordinatedPhylogeneticConservatism),
                                                                 "Coordinated_phylogenetic_niche_conservatism" = mean(coordinatedNichePhylogeneticConservatism),
                                                                 "Total_environmental_coordination" = mean(totalEnvironmentalCoordination),
                                                                 "Pure_environmental_coordination" = mean(pureEnvironmentalCoordination),
                                                                 "Residual_coordination" = mean(residualCoordination),

                                                                 "Pvalue_Pure_coordinated_phylogenetic_conservatism" = pureCoordinatedPhylogeneticConservatism_pval,
                                                                 "Pvalue_Coordinated_phylogenetic_niche_conservatism" = coordinatedNichePhylogeneticConservatism_pval,
                                                                 "Pvalue_Total_environmental_coordination" = totalEnvironmentalCoordination_pval,
                                                                 "Pvalue_Pure_environmental_coordination" = pureEnvironmentalCoordination_pval,
                                                                 "Pvalue_Residual_coordination" = residualCoordination_pval
        )

        covariancePartitionResults$covariancePartitionDistributions[["pureCoordinatedPhylogeneticConservatism"]] <- pureCoordinatedPhylogeneticConservatism
        covariancePartitionResults$covariancePartitionDistributions[["coordinatedPhylogeneticNicheConservatism"]] <- coordinatedNichePhylogeneticConservatism
        covariancePartitionResults$covariancePartitionDistributions[["totalEnvironmentalCoordination"]] <- totalEnvironmentalCoordination
        covariancePartitionResults$covariancePartitionDistributions[["pureEnvironmentalCoordination"]] <- pureEnvironmentalCoordination
        covariancePartitionResults$covariancePartitionDistributions[["residualCoordination"]] <- residualCoordination

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

  if(save){
    if(!is.null(environmental.variables)){
      save(list = paste0("traitsCovariancePartitionResults"), file = paste0(outputs.dir, "/models_outputs/partialTraitsCovariancePartitionResults.RData"))
      print(paste0(outputs.dir, "/models_outputs/partialTraitsCovariancePartitionResults.RData"))
    } else{
      save(list = paste0("traitsVariancePartitionResults"), file = paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults.RData"))
      print(paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults.RData"))
    }
  }
  return(traitsCovariancePartitionResults)
}
