#' Variance partition including (or not) environment
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
computeVariancePartition <- function(traits, environmental.variables = NULL, dataset, phylogeny, model.specifications = NULL, force.run = T, save = T) {

  # results object
  traitsVariancePartitionResults <- list()
  traitsVariancePartitionResults$varianceResults <- data.frame()
  traitsVariancePartitionResults$models.diagnostics <- data.frame()
  traitsVariancePartitionResults$individual.models.results <- list()

  # prepare models structure
  uni_mdls.str <- data.frame(trait = traits)
  uni_mdls.str$type <- paste0("uni_", uni_mdls.str$trait)
  uni_mdls.str$pred_var <- ""
  uni_mdls.str$fix.frml <- paste0(uni_mdls.str$trait, " ~ 1")
  uni_mdls.str$ran.frml <- "~ animal"

  # lad previous results, if exist

  if(!is.null(results.file)){
    results.file <- paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults", environmental.variables, ".RData")
  } else {
    results.file <- paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults.RData")
  }

  if (file.exists(results.file) && isFALSE(force.run)) {
    print("loanding previous results")
    load(file = results.file)
  }

  # run models and extract results
  for (model in uni_mdls.str$type) {

    if (!model %in% names(traitsVariancePartitionResults$individual.models.results) | force.run) {

      print(paste0("Running variance calculation: ", model))

      model.descr <- uni_mdls.str %>%
        dplyr::filter(type == model)

      trait <- model.descr$trait

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait, environmental.variables))

      if (is.null(model.specifications)) {
        print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
        model.specifications <- defineModelsSpecifications()
      }


      ### Variance partition without accounting for environment ####

      ## Phylogenetic model

      fix.frml <- paste0(trait, "~ 1",  collapse = " ")

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

      totalNonPhylogenetic <- totalNonPhyloVar / totalVar

      # results
      variancePartitionResults <- list()
      variancePartitionResults$variancePartition <- data.frame("Trait" = trait,
                                                               "N" = length(modellingData$dta$animal),
                                                               "Total_phylogenetic_conservatism" = mean(totalPhylogeneticConservatism),
                                                               "Total_non_phylogenetic" = mean(totalNonPhylogenetic)
      )

      variancePartitionResults$variancePartitionDistributions <- list("totalPhylogeneticConservatism" = totalPhylogeneticConservatism,
                                                                      "totalNonPhylogenetic" = totalNonPhylogenetic)
      variancePartitionResults$modelPhylo <- mdlPhylo
      variancePartitionResults$model.diagnostics <- model.diagnostics

      if(is.null(environmental.variables)){


      # add to all traits results
      traitsVariancePartitionResults$varianceResults <- rbind(traitsVariancePartitionResults$varianceResults,
                                                                          variancePartitionResults$variancePartition)

      traitsVariancePartitionResults$models.diagnostics <- rbind(traitsVariancePartitionResults$models.diagnostics,
                                                                 variancePartitionResults$modelPhylo.diagnostics)
      traitsVariancePartitionResults$individual.models.results[[model]] <- variancePartitionResults
}

      ### Variance partition including the environment ####

      if(!is.null(environmental.variables)){

        for(predictor in environmental.variables){
          fix.frml <- paste0(fix.frml, " + ", predictor, collapse = " ")
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
                                                            "Environmental_variables" = paste0(environmental.variables, collapse = ", "),
                                                            "Pure_phylogenetic_conservatism" = mean(purePhylogeneticConservatism),
                                                            "Phylogenetic_niche_conservatism" = mean(phylogeneticNicheConservatism),
                                                            "Total_environmental" = mean(totalEnvironmental),
                                                            "Pure_environmental" = mean(pureEnvironmental),
                                                            "Residual" = mean(residual)
        )

        variancePartitionResults$variancePartitionDistributions[["Pure_phylogenetic_conservatism"]] <- purePhylogeneticConservatism
        variancePartitionResults$variancePartitionDistributions[["Phylogenetic_niche_conservatism"]] <- phylogeneticNicheConservatism
        variancePartitionResults$variancePartitionDistributions[["total_environmental"]] <- totalEnvironmental
        variancePartitionResults$variancePartitionDistributions[["pure_environmental"]] <- pureEnvironmental
        variancePartitionResults$variancePartitionDistributions[["residual"]] <- residual

        variancePartitionResults$modelPhyloEnv <- mdlPhyloEnv
        variancePartitionResults$model.diagnostics <- rbind(variancePartitionResults$model.diagnostics, model.diagnostics)

        # add to all traits results
        traitsVariancePartitionResults$varianceResults <- rbind(traitsVariancePartitionResults$varianceResults,
                                                                            variancePartitionResults$variancePartition)

        traitsVariancePartitionResults$models.diagnostics <- rbind(traitsVariancePartitionResults$models.diagnostics,
                                                                   variancePartitionResults$model.diagnostics)
        traitsVariancePartitionResults$individual.models.results[[model]] <- variancePartitionResults

      }

    }  # end evaluation if model already exists in results

  } # end bucle for all traits

  print("Model structure used:")
  print(uni_mdls.str)
  print("Phylogenetic signal results:")
  print(traitsVariancePartitionResults$varianceResults)

  # save results

  if(save){
      save(list = paste0("traitsVariancePartitionResults"), file = results.file)
      print(results.file)
  }
  return(traitsVariancePartitionResults)
}
