#' Variance partition including (or not) environment
#'
#' @param trait (character) Name of the trait or list of traits. It  must be contained in the dataset.
#' @param environmentalVariables (character) Names of the environmental variables They must be contained in dataset.
#' @param dataset (data frame) Dataset containing the trait of interest and a column named "animal" describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param model.specifications (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @param forceRun (logical) If false, models already run are not runned again.
#'
#' @return
#' @export
#'
#' @examples
computeVariancePartition <- function(traits, environmentalVariables = NULL, dataset, phylogeny, model.specifications = NULL, forceRun = F) {

  # results object
  traitsVariancePartitionResults <- list()
  traitsVariancePartitionResults$phylogenetic.signal.results <- data.frame()
  traitsVariancePartitionResults$models.diagnostics <- data.frame()
  traitsVariancePartitionResults$individual.models.results <- list()

  # prepare models structure
  uni_mdls.str <- data.frame(resp_var = traits)
  uni_mdls.str$type <- paste0("uni_", uni_mdls.str$resp_var)
  uni_mdls.str$n_respVars <- 1
  uni_mdls.str$pred_var <- ""
  uni_mdls.str$fix.frml <- paste0(uni_mdls.str$resp_var, " ~ 1")
  uni_mdls.str$ran.frml <- "~ animal"

  # lad previous results, if exist
  if (file.exists(paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults.RData")) && isFALSE(FORCERUN)) {
    print("loanding previous results")
    load(file = paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults.RData"))
  }

  # run models and extract results
  for (model in uni_mdls.str$type) {

    if (!model %in% names(traitsVariancePartitionResults$individual.models.results) | FORCERUN) {

      print(paste0("Running phylo. signal model: ", model))

      model.descr <- uni_mdls.str %>%
<<<<<<< HEAD
        dplyr::filter(type == model)
=======
                     dplyr::filter(type == model)
>>>>>>> bf0db1d746dc56f588e94da480ef4be98f4a6cc3

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait, environmentalVariables))

      if (is.null(model.specifications)) {
        print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
        model.specifications <- defineModelsSpecifications()
      }


      ### Variance partition without accounting for environment ####

      ## Phylogenetic model

      fix.frml <- paste0(trait, "~ 1 ")

      mdlPhylo <- MCMCglmm::MCMCglmm(stats::as.formula(fix.frml),
                                     random = ~animal, family = "gaussian", prior = model.specifications$uniresponse_prior,
                                     data = modellingData$dta, pedigree = modellingData$phylo,
                                     nitt = model.specifications$number_interations, burnin = model.specifications$burning_iterations,
                                     thin = model.specifications$thinning_iterations,
                                     verbose = F)

      mdlPhylo$name <- fix.frml

      model.diagnostics <- diagnoseModels(model = mdlPhylo)


      ### Variance partition calculation

      # total variance
      totalVar <-  mdlPhylo$VCV[, "animal"] +  mdlPhylo$VCV[, "units"]

      # total phylogenetic variance
      totalPhyloVar <- mdlPhylo$VCV[, "animal"]

      # total residual variance
      totalNonPhyloVar <- mdlPhylo$VCV[, "units"]


      ## results

      # total phylogenetic conservatism
      totalPhylogeneticConservatism <-  totalPhyloVar / totalVar

      # results
      variancePartitionResults <- list()
      variancePartitionResults$variancePartition <- data.frame("Trait" = trait,
                                                               "N" = length(modellingData$dta$animal),
                                                               "Total_phylogenetic_conservatism" = mean(totalPhylogeneticConservatism)
      )
      variancePartitionResults$variancePartitionDistributions <- list("totalPhylogeneticConservatism" = totalPhylogeneticConservatism)
      variancePartitionResults$modelPhylo <- mdlPhylo
      variancePartitionResults$modelPhylo.diagnostics <- model.diagnostics


      ### Variance partition including the environment ####

      if(!is.null(environmentalVariables)){

        for(predictor in environmentalVariables){
          fix.frml <- paste0(fix.frml, " + ", predictor)
        }

        mdlPhyloEnv <- MCMCglmm::MCMCglmm(stats::as.formula(fix.frml),
                                          random = ~animal, family = "gaussian", prior = model.specifications$uniresponse_prior,
                                          data = modellingData$dta, pedigree = modellingData$phylo,
                                          nitt = model.specifications$number_interations, burnin = model.specifications$burning_iterations,
                                          thin = model.specifications$thinning_iterations,
                                          verbose = F)

        mdlPhyloEnv$name <- fix.frml

        model.diagnostics <- diagnoseModels(model = mdlPhyloEnv)

        ### Variance partition calculation

        # non-environmental phylogenetic variance (pure phylogenetic variance)
        purePhyloVar <- mdlPhyloEnv$VCV[, "animal"]

        # avoid negative variances when purePhylo variance is a little bit higher than totalPhylo variance
        if(mean(purePhyloVar) > mean(totalPhyloVar)){
          purePhyloVar <- totalPhyloVar
        }

        # Residual variance
        pureResidualVar <- mdlPhyloEnv$VCV[, "units"]

        ## results

        # pure phylogenetic conservatism
        purePhylogeneticConservatism <-  purePhyloVar / totalVar

        # phylogenetic niche conservatism (environment x phylogeny)
        phylogeneticNicheConservatism <- (totalPhyloVar - purePhyloVar) / totalVar

        # pure environmental
        pureEnvironmental <- (totalNonPhyloVar - pureResidualVar) / totalVar

        # total environmental
        totalEnvironmental <- pureEnvironmental + phylogeneticNicheConservatism

        # residual
        residual <- pureResidualVar / totalVar


        # results

        variancePartitionResults$variancePartition <- cbind(variancePartitionResults$variancePartition,
                                                            "Environmental_variables" = paste0(environmentalVariables, collapse = ", "),
                                                            "Pure_phylogenetic_conservatism" = mean(purePhyloVar),
                                                            "Phylogenetic_niche_conservatism" = mean(phylogeneticNicheConservatism),
                                                            "total_environmental" = mean(totalEnvironmental),
                                                            "pure_environmental" = mean(pureEnvironmental),
                                                            "residual" = mean(residual)
<<<<<<< HEAD
        )

        variancePartitionResults$variancePartitionDistributions[Pure_phylogenetic_conservatism] <- purePhyloVar
        variancePartitionResults$variancePartitionDistributions[Phylogenetic_niche_conservatism] <- phylogeneticNicheConservatism
        variancePartitionResults$variancePartitionDistributions[total_environmental] <- totalEnvironmental
        variancePartitionResults$variancePartitionDistributions[pure_environmental] <- pureEnvironmental
        variancePartitionResults$variancePartitionDistributions[residual] <- residual

        variancePartitionResults$mdlPhyloEnv <- mdlPhyloEnv
        variancePartitionResults$mdlPhyloEnv.diagnostics <- model.diagnostics

        # add to all traits results
        traitsVariancePartitionResults$phylogenetic.signal.results <- rbind(traitsVariancePartitionResults$phylogenetic.signal.results,
                                                                            variancePartitionResults$phyloSignal)

        traitsVariancePartitionResults$models.diagnostics <- rbind(traitsVariancePartitionResults$models.diagnostics,
                                                                   variancePartitionResults$model.diagnostics)
        traitsVariancePartitionResults$individual.models.results[[model]] <- variancePartitionResults

      }
=======
                                                            )

    variancePartitionResults$variancePartitionDistributions[Pure_phylogenetic_conservatism] <- purePhyloVar
    variancePartitionResults$variancePartitionDistributions[Phylogenetic_niche_conservatism] <- phylogeneticNicheConservatism
    variancePartitionResults$variancePartitionDistributions[total_environmental] <- totalEnvironmental
    variancePartitionResults$variancePartitionDistributions[pure_environmental] <- pureEnvironmental
    variancePartitionResults$variancePartitionDistributions[residual] <- residual

    variancePartitionResults$mdlPhyloEnv <- mdlPhyloEnv
    variancePartitionResults$mdlPhyloEnv.diagnostics <- model.diagnostics

    # add to all traits results
    traitsVariancePartitionResults$phylogenetic.signal.results <- rbind(traitsVariancePartitionResults$phylogenetic.signal.results,
                                                                   variancePartitionResults$phyloSignal)

    traitsVariancePartitionResults$models.diagnostics <- rbind(traitsVariancePartitionResults$models.diagnostics,
                                                          variancePartitionResults$model.diagnostics)
    traitsVariancePartitionResults$individual.models.results[[model]] <- variancePartitionResults

  }
>>>>>>> bf0db1d746dc56f588e94da480ef4be98f4a6cc3

    }  # end evaluation if model already exists in results

  } # end bucle for all traits

  return(traitsVariancePartitionResults)
}
