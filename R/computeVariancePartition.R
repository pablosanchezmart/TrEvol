#' Phylogenetic variance partition including (or not) one or more environmental variables
#'
#'Compute total variance and separate it into its phylogenetic and non-phylogenetic components.
#' If environmental variable is included, four variance is separated into four components: non-attributed phylogenetic variance, environmental
#'phylogenetic variance, non-phylogenetic (labile) environmental variance and residual variance. Credible intervals for each estimate are reported,
#'from which statistical significance can be assessed (significant when CI does not contain zero). Significance based on Credible Intervals is displayed for each estimate. P-value related
#'to statistical significance are also provided as an alternative assessment.
#'This is a version of the computeVariancevariancePartition function that allows to include more than one environmental factor at the same time.
#'
#'The function can automaticallys save results. For that, first run initializeTrEvol to create the folders and subfolders structure.
#'
#' @param traits (*character*). Name of the trait or list of traits to calculate variances. They must be in the dataset.
#' @param environmental_variables (*character*). Names of the environmental variables They must be contained in dataset.
#' @param dataset (*data frame*). Data frame containing the trait of interest and a column describing terminal appearing as tip labels of the phylogeny.
#' @param terminal_taxa (*character*). Terminal taxon as named in the dataset (e.g., species).
#' @param phylogeny (*phylo*) Phylogeny with tip labels contained in the terminal taxon column of the dataset.
#' @param model_specifications (*list*). Mcmcglmm models specifications as specified by the defineModelsSpecification of this package. If not defined, the function internally uses the default in defineModelsSpecification function.
#' @param force_run (*logical*) If false, models already run and saved in outputs folder are not run again.
#' @param save (*logical*) If false, results are not saved in the outputs folder.
#' @param verbose (*logical*) If false, function runs quietly. Otherwise, iterations and diagnostics are shown while running.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate example data
#' simulated_traits.data <- simulateDataSet()
#'
#' # Compute variance structure for simulated traits using default parameters
#' variance_results <- computeVariancePartition(
#' traits = c("phylo_G1_trait1", "phylo_G1_trait2"),
#' environmental_variable = c("phylo_G1_env", "phylo_G2_env")
#' data = simulated_traits.data$data
#' phylogeny = simulated_traits.data$phylogeny
#' )
#' }
#'
computeVariancePartition <- function(traits = NULL,
                                     environmental_variables = NULL,
                                     terminal_taxa = NULL,
                                     dataset = NULL,
                                     phylogeny = NULL,
                                     model_specifications = NULL,
                                     force_run = F,
                                     save = F,
                                     verbose = F) {



  # Arguments
  if(is.null(traits)){
    stop("Specify traits argument")
  }

  if(is.null(dataset)){
    stop("Specify dataset argument")
  }

  if(is.null(terminal_taxa)){
    stop("Specify terminal_taxa argument")
  }

  if(is.null(phylogeny)){
    stop("Specify phylogeny argument")
  }


  # Check whether initializeTrevol has been run if user wants to automatically save results
  if(isTRUE(save)){
    if(isFALSE(exists("outputs.dir"))){
      stop("If you set save = T you first need to run initializeTrEvol to create the folder and subfolder structure where results are saved.")
    }
  }

  # Make sure that tip nodes are deleted as they can give some problems
  if(!is.null(phylogeny$tip.nodes)){
    phylogeny$tip.nodes <- NULL
  }


  # Name a terminal taxon column as animal for MCMCglmm

  dataset$animal <- dataset[, terminal_taxa]

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

  if(!is.null(environmental_variables)){
    results.file <- paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults_", paste0(environmental_variables, collapse = "_"), ".RData")
  } else {
    results.file <- paste0(outputs.dir, "/models_outputs/traitsVariancePartitionResults.RData")
  }

  if (file.exists(results.file) && isFALSE(force_run)) {
    message("loanding previous results")
    load(file = results.file)
  }

  # run models and extract results
  for (model in uni_mdls.str$type) {

    if (!model %in% names(traitsVariancePartitionResults$individual.models.results) | force_run) {

      message(paste0("Running variance calculation: ", model))

      model.descr <- uni_mdls.str %>%
        dplyr::filter(type == model)

      trait <- model.descr$trait

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait, environmental_variables))

      if (is.null(model_specifications)) {
        message("Using default model specificatios. Use defineModelsSpecifications() output on model_specifications argument to set them manually.")
        model_specifications <- defineModelsSpecifications()
      }


      ### Variance partition without accounting for environment ####

      ## Phylogenetic model

      fix.frml <- paste0(trait, "~ 1",  collapse = " ")

      mdlPhylo <- MCMCglmm::MCMCglmm(stats::as.formula(fix.frml),
                                     random = ~animal, family = "gaussian",
                                     prior = model_specifications$uniresponse_prior,
                                     data = modellingData$dta, pedigree = modellingData$phylo,
                                     nitt = model_specifications$number_iterations,
                                     burnin = model_specifications$burning_iterations,
                                     thin = model_specifications$thinning_iterations,
                                     verbose = verbose)

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
      total_phylogenetic_variance_t1 <-  totalPhyloVar / totalVar

      # Significance
      total_phylogenetic_variance_t1_lCI <- as.numeric(RChronoModel::CredibleInterval(total_phylogenetic_variance_t1)[2])
      total_phylogenetic_variance_t1_hCI <- as.numeric(RChronoModel::CredibleInterval(total_phylogenetic_variance_t1)[3])

      total_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_variance_t1))) # Pvalue


      total_non_phylogenetic_variance_t1 <- totalNonPhyloVar / totalVar

      # Significance
      total_non_phylogenetic_variance_t1_lCI <- as.numeric(RChronoModel::CredibleInterval(total_non_phylogenetic_variance_t1)[2])
      total_non_phylogenetic_variance_t1_hCI <- as.numeric(RChronoModel::CredibleInterval(total_non_phylogenetic_variance_t1)[3])


      total_non_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_variance_t1))) # Pvalue

      # results
      variancePartitionResults <- list()
      variancePartitionResults$variancePartition <- data.frame("trait" = trait,
                                                               "number_observations" = length(modellingData$dta$animal),

                                                               "phylogenetic_variance" = mean(total_phylogenetic_variance_t1),
                                                               "non_phylogenetic_variance" = mean(total_non_phylogenetic_variance_t1),

                                                               "p_value_phylogenetic_variance" = total_phylogenetic_variance_t1_pvalue,
                                                               "low_CI_phylogenetic_variance" = total_phylogenetic_variance_t1_lCI,
                                                               "high_CI_phylogenetic_variance" = total_phylogenetic_variance_t1_hCI,
                                                               "CI_significance_phylogenetic_variance" = ifelse(total_phylogenetic_variance_t1_lCI < 0 && total_phylogenetic_variance_t1_hCI > 0, "no", "yes")
      )

      variancePartitionResults$variancePartitionDistributions[[paste0("phylogenetic_variance_", trait)]] <- total_phylogenetic_variance_t1
      variancePartitionResults$variancePartitionDistributions[[paste0("non_phylogenetic_variance_", trait)]] <- total_non_phylogenetic_variance_t1

      variancePartitionResults$modelPhylo <- mdlPhylo
      variancePartitionResults$model.diagnostics <- model.diagnostics

      if(is.null(environmental_variables)){


      # add to all traits results
      traitsVariancePartitionResults$varianceResults <- rbind(traitsVariancePartitionResults$varianceResults,
                                                                          variancePartitionResults$variancePartition)

      traitsVariancePartitionResults$models.diagnostics <- rbind(traitsVariancePartitionResults$models.diagnostics,
                                                                 variancePartitionResults$modelPhylo.diagnostics)
      traitsVariancePartitionResults$individual.models.results[[model]] <- variancePartitionResults
}

      ### Variance partition including the environment ####

      if(!is.null(environmental_variables)){

        for(predictor in environmental_variables){
          fix.frml <- paste0(fix.frml, " + ", predictor, collapse = " ")
        }

        mdlPhyloEnv <- MCMCglmm::MCMCglmm(stats::as.formula(fix.frml),
                                          random = ~animal, family = "gaussian",
                                          prior = model_specifications$uniresponse_prior,
                                          data = modellingData$dta, pedigree = modellingData$phylo,
                                          nitt = model_specifications$number_iterations,
                                          burnin = model_specifications$burning_iterations,
                                          thin = model_specifications$thinning_iterations,
                                          verbose = verbose)

        mdlPhyloEnv$name <- fix.frml

        model.diagnostics <- diagnoseModels(model = mdlPhyloEnv)

        ### Variance partition calculation

        # total variance

        totalVar2 <- mdlPhyloEnv$VCV[, "animal"] + mdlPhyloEnv$VCV[, "units"]

        # non-environmental phylogenetic variance (pure phylogenetic variance)
        purePhyloVar <- mdlPhyloEnv$VCV[, "animal"]

        # avoid negative variances when purePhylo variance is a little bit higher than totalPhylo variance
        if(mean(purePhyloVar) > mean(totalPhyloVar)){
          purePhyloVar <- totalPhyloVar
        }

        # Residual variance
        pureResidualVar <- mdlPhyloEnv$VCV[, "units"]

        ## results

        # non-attributed phylogenetic conservatism

        non_attributed_phylogenetic_variance_t1 <-  purePhyloVar / totalVar2

        # Significance
        non_attributed_phylogenetic_variance_t1_lCI <- as.numeric(RChronoModel::CredibleInterval(non_attributed_phylogenetic_variance_t1)[2])
        non_attributed_phylogenetic_variance_t1_hCI <- as.numeric(RChronoModel::CredibleInterval(non_attributed_phylogenetic_variance_t1)[3])

        non_attributed_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(non_attributed_phylogenetic_variance_t1))) # Pvalue



        # phylogenetic niche conservatism (environment x phylogeny)
        environmental_phylogenetic_variance_t1 <- (totalPhyloVar - purePhyloVar) / totalVar2

        # Significance
        environmental_phylogenetic_variance_t1_lCI <- as.numeric(RChronoModel::CredibleInterval(environmental_phylogenetic_variance_t1)[2])
        environmental_phylogenetic_variance_t1_hCI <- as.numeric(RChronoModel::CredibleInterval(environmental_phylogenetic_variance_t1)[3])

        environmental_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_variance_t1))) # Pvalue


        # labile environmental
        labile_environmental_variance_t1 <- (totalNonPhyloVar - pureResidualVar) / totalVar2

        # Significance
        labile_environmental_variance_t1_lCI <- as.numeric(RChronoModel::CredibleInterval(labile_environmental_variance_t1)[2])
        labile_environmental_variance_t1_hCI <- as.numeric(RChronoModel::CredibleInterval(labile_environmental_variance_t1)[3])

        labile_environmental_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(labile_environmental_variance_t1))) # Pvalue


        # total environmental
        totalEnvironmental <- labile_environmental_variance_t1 + environmental_phylogenetic_variance_t1

        # residual
        residual_variance_t1 <- pureResidualVar / totalVar2

        # Significance
        residual_variance_t1_lCI <- as.numeric(RChronoModel::CredibleInterval(residual_variance_t1)[2])
        residual_variance_t1_hCI <- as.numeric(RChronoModel::CredibleInterval(residual_variance_t1)[3])

        residual_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(residual_variance_t1))) # Pvalue

        # results

        variancePartitionResults$variancePartition <- cbind(variancePartitionResults$variancePartition,
                                                            "environmental_variable" = paste0(environmental_variables, collapse = ", "),
                                                            "non_attributed_phylogenetic_variance" = mean(non_attributed_phylogenetic_variance_t1),
                                                            "environmental_phylogenetic_variance" = mean(environmental_phylogenetic_variance_t1),
                                                            # "Total_environmental" = mean(totalEnvironmental),
                                                            "labile_environmental_variance" = mean(labile_environmental_variance_t1),
                                                            "residual_variance" = mean(residual_variance_t1),

                                                            "p_value_non_attributed_phylogenetic_variance" = non_attributed_phylogenetic_variance_t1_pvalue,
                                                            "low_CI_non_attributed_phylogenetic_variance" = non_attributed_phylogenetic_variance_t1_lCI,
                                                            "high_CI_non_attributed_phylogenetic_variance" = non_attributed_phylogenetic_variance_t1_hCI,
                                                            "CI_significance_non_attributed_phylogenetic_variance" = ifelse(non_attributed_phylogenetic_variance_t1_lCI < 0 && non_attributed_phylogenetic_variance_t1_hCI > 0, "no", "yes"),

                                                            "p_value_environmental_phylogenetic_variance" = environmental_phylogenetic_variance_t1_pvalue,
                                                            "low_CI_environmental_phylogenetic_variance" = environmental_phylogenetic_variance_t1_lCI,
                                                            "high_CI_environmental_phylogenetic_variance" = environmental_phylogenetic_variance_t1_hCI,
                                                            "CI_significance_environmental_phylogenetic_variance" = ifelse(environmental_phylogenetic_variance_t1_lCI < 0 && environmental_phylogenetic_variance_t1_hCI > 0, "no", "yes"),

                                                            "p_value_labile_environmental_variance" = labile_environmental_variance_t1_pvalue,
                                                            "low_CI_labile_environmental_variance" = labile_environmental_variance_t1_lCI,
                                                            "high_CI_labile_environmental_variance" = labile_environmental_variance_t1_hCI,
                                                            "CI_significance_labile_environmental_variance" = ifelse(labile_environmental_variance_t1_lCI < 0 && labile_environmental_variance_t1_hCI > 0, "no", "yes")
        )

        variancePartitionResults$variancePartitionDistributions[[paste0("non_attributed_phylogenetic_variance", trait)]] <- non_attributed_phylogenetic_variance_t1
        variancePartitionResults$variancePartitionDistributions[[paste0("environmental_phylogenetic_variance", trait)]] <- environmental_phylogenetic_variance_t1
        variancePartitionResults$variancePartitionDistributions[[paste0("labile_environmental_variance", trait)]] <- labile_environmental_variance_t1
        variancePartitionResults$variancePartitionDistributions[[paste0("residual_variance", trait)]] <- residual_variance_t1


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

  # message("Model structure used:")
  # message(uni_mdls.str)
  # message("Phylogenetic signal results:")
  # message(traitsVariancePartitionResults$varianceResults)

  # save results

  if(save){
    if(!is.null(environmental_variables)){
      assign(paste0("traitsVariancePartitionResults_", environmental_variables, collapse = "_"), traitsVariancePartitionResults)
      save(list = paste0("traitsVariancePartitionResults_", environmental_variables, collapse = "_"), file = results.file)
      message(results.file)
    } else{
      save(list = paste0("traitsVariancePartitionResults"), file = results.file)
      message(results.file)
    }
  }
  return(traitsVariancePartitionResults)
}
