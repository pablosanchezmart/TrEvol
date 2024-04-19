#' Phylogenetic covariance partition including (or not) one or more environmental variables
#'
#'Compute total covariance and separate it into its phylogenetic and non-phylogenetic components.
#' If environmental variable is included, four covariance is separated into four components: non-attributed phylogenetic covariance, environmental
#'phylogenetic covariance, non-phylogenetic (labile) environmental covariance and residual covariance. Credible intervals for each estimate are reported,
#'from which statistical significance can be assessed (significant when CI does not contain zero). Significance based on Credible Intervals is displayed for each estimate. P-value related
#'to statistical significance are also provided as an alternative assessment.
#'This is a version of the computeVarianceCovariancePartition function that allows to include more than one environmental factor at the same time.
#'
#'The function can automaticallys save results. For that, first run initializeTrEvol to create the folders and subfolders structure.
#'
#' @param traits (*character*). Name of the trait or list of traits to calculate covariances. They must be in the dataset.
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
#' # Compute covariance structure for simulated traits using default parameters
#' covariance_results <- computeCovariancePartition(
#' traits = c("phylo_G1_trait1", "phylo_G1_trait2"),
#' environmental_variable = c("phylo_G1_env", "phylo_G2_env")
#' data = simulated_traits.data$data
#' phylogeny = simulated_traits.data$phylogeny
#' )
#' }

computeCovariancePartition <- function(traits = NULL,
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
  traitsCVPartitionResults <- list()
  traitsCVPartitionResults$covariancePartition <- data.frame()
  traitsCVPartitionResults$model.diagnostics <- data.frame()
  traitsCVPartitionResults$individual.models.results <- list()


  ### MODELS STRUCTURE --------------------------------------------------------- ####

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

  # load previous results, if exist

  if(!is.null(environmental_variables)){
    results.file <- paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults", paste0(environmental_variables, collapse = "_"), ".RData")
  } else {
    results.file <- paste0(outputs.dir, "/models_outputs/traitsCovariancePartitionResults.RData")
  }

  if (file.exists(results.file) && isFALSE(force_run)) {
    message("loanding previous results")
    load(file = results.file)
  }


  ### RUN MODELS AND EXTRACT RESULTS ------------------------------------------- ####

  for (model in multi_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(traitsCVPartitionResults$individual.models.results) | force_run) {

      message(paste0("Running covariance calculation: ", model))
      model.descr <- multi_mdls.str %>%
        dplyr::filter(type == model)

      trait1 <- as.character(model.descr$trait1)
      trait2 <- as.character(model.descr$trait2)

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait1, trait2, environmental_variables))

      # formula
      fix.frml <- paste0("cbind(", trait1, ", ", trait2, ") ~ trait-1")

      if (is.null(model_specifications)) {
        message("Using default model specificatios. Use defineModelsSpecifications() output on model_specifications argument to set them manually.")
        model_specifications <- defineModelsSpecifications()
      }

      ### Covariance partition without accounting for environment ####

      ## Phylogenetic model

      mdlPhylo <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                     random = ~ us(trait):animal, rcov = ~us(trait):units,
                                     data= modellingData$dta, pedigree = modellingData$phylo,
                                     family = c("gaussian", "gaussian"),
                                     prior = model_specifications$biresponse_prior,
                                     nitt = model_specifications$number_iterations,
                                     burnin = model_specifications$burning_iterations,
                                     thin = model_specifications$thinning_iterations,
                                     verbose = verbose)

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

      total_covariance  <- totalCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ) )

      # Significance
      total_covariance_lCI <- as.numeric(RChronoModel::CredibleInterval(total_covariance)[2])
      total_covariance_hCI <- as.numeric(RChronoModel::CredibleInterval(total_covariance)[3])

      total_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_covariance))) # Pvalue

      # total coordinated phylogenetic conservatism

      total_phylogenetic_covariance <- totalPhyloCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

      # Significance
      total_phylogenetic_covariance_lCI <- as.numeric(RChronoModel::CredibleInterval(total_phylogenetic_covariance)[2])
      total_phylogenetic_covariance_hCI <- as.numeric(RChronoModel::CredibleInterval(total_phylogenetic_covariance)[3])

      total_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_covariance))) # Pvalue

      # total non-phylogenetic coordination

      total_non_phylogenetic_covariance  <- totalResidualCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

      # Significance
      total_non_phylogenetic_covariance_lCI <- as.numeric(RChronoModel::CredibleInterval(total_non_phylogenetic_covariance)[2])
      total_non_phylogenetic_covariance_hCI <- as.numeric(RChronoModel::CredibleInterval(total_non_phylogenetic_covariance)[3])

      total_non_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_covariance))) # Pvalue

      CVPartitionResults <- list()

      CVPartitionResults$covariancePartition <-  data.frame("trait_1" = trait1,
                                                                    "trait_2" = trait2,
                                                                    "number_observations" = length(modellingData$dta$animal),

                                                                    "total_correlation" = mean(total_covariance),
                                                                    "phylogenetic_correlation" = mean(total_phylogenetic_covariance),
                                                                    "non_phylogenetic_correlation" = mean(total_non_phylogenetic_covariance),

                                                                    "p_value_total_correlation" = total_covariance_pvalue,
                                                                    "low_CI_total_correlation" = total_covariance_lCI,
                                                                    "high_CI_total_correlation" = total_covariance_hCI,
                                                                    "CI_significance_total_correlation" = ifelse(total_covariance_lCI < 0 && total_covariance_hCI > 0, "no", "yes"),


                                                                    "p_value_phylogenetic_correlation" = total_phylogenetic_covariance_pvalue,
                                                                    "low_CI_phylogenetic_correlation" = total_phylogenetic_covariance_lCI,
                                                                    "high_CI_phylogenetic_correlation" = total_phylogenetic_covariance_hCI,
                                                                    "CI_significance_phylogenetic_correlation" = ifelse(total_phylogenetic_covariance_lCI < 0 && total_phylogenetic_covariance_hCI > 0, "no", "yes"),

                                                                    "p_value_non_phylogenetic_correlation" = total_non_phylogenetic_covariance_pvalue,
                                                                    "low_CI_non_phylogenetic_correlation" = total_non_phylogenetic_covariance_lCI,
                                                                    "high_CI_non_phylogenetic_correlation" = total_non_phylogenetic_covariance_hCI,
                                                                    "CI_significance_non_phylogenetic_correlation" = ifelse(total_non_phylogenetic_covariance_lCI < 0 && total_non_phylogenetic_covariance_hCI > 0, "no", "yes")
      )

      CVPartitionResults$covariancePartitionDistributions <- list("total_correlation" = total_covariance,
                                                                   "phylogenetic_correlation" = total_phylogenetic_covariance,
                                                                   "non_phylogenetic_correlation" = total_non_phylogenetic_covariance)

      CVPartitionResults$modelPhylo <- mdlPhylo
      CVPartitionResults$model.diagnostics <- model.diagnostics


      if(is.null(environmental_variables)){
        # add to all traits results
        traitsCVPartitionResults$covariancePartition <- rbind(traitsCVPartitionResults$covariancePartition,
                                                                    CVPartitionResults$covariancePartition)

        traitsCVPartitionResults$model.diagnostics <- rbind(traitsCVPartitionResults$model.diagnostics,
                                                                     CVPartitionResults$model.diagnostics)
        traitsCVPartitionResults$individual.models.results[[model]] <- CVPartitionResults
      }


      ### Variance partition including the environment ####

      if(!is.null(environmental_variables)){

        for(predictor in environmental_variables){
          fix.frml <- paste0(fix.frml, " + trait:", predictor)
        }

        mdlPhyloEnv <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                          random = ~ us(trait):animal, rcov = ~us(trait):units,
                                          data= modellingData$dta, pedigree = modellingData$phylo,
                                          family = c("gaussian", "gaussian"),
                                          prior = model_specifications$biresponse_prior,
                                          nitt = model_specifications$number_iterations,
                                          burnin = model_specifications$burning_iterations,
                                          thin = model_specifications$thinning_iterations,
                                          verbose = verbose)

        mdlPhyloEnv$name <- fix.frml

        model.diagnostics <- diagnoseModels(model = mdlPhyloEnv)


        ### Covariance partition calculation

        # non-attributed phylogenetic covariance
        purePhyloCV <-  mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")]

        # avoid negative variances when purePhylo covariance is a little bit higher than totalPhylo variance
        # if(mean(purePhyloCV) > mean(totalPhyloCV)){
        #   purePhyloCV <- totalPhyloCV
        # }

        # pure residual covariance
        pureResidualCV <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]


        # total coordinated phylogenetic conservatism

        non_attributed_phylogenetic_covariance <- purePhyloCV / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

        # Significance
        non_attributed_phylogenetic_covariance_lCI <- as.numeric(RChronoModel::CredibleInterval(non_attributed_phylogenetic_covariance)[2])
        non_attributed_phylogenetic_covariance_hCI <- as.numeric(RChronoModel::CredibleInterval(non_attributed_phylogenetic_covariance)[3])

        non_attributed_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(non_attributed_phylogenetic_covariance))) # Pvalue


        # total coordinated phylogenetic niche conservatism

        environmental_phylogenetic_covariance <- (totalPhyloCV - purePhyloCV) / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

        # Significance
        environmental_phylogenetic_covariance_lCI <- as.numeric(RChronoModel::CredibleInterval(environmental_phylogenetic_covariance)[2])
        environmental_phylogenetic_covariance_hCI <- as.numeric(RChronoModel::CredibleInterval(environmental_phylogenetic_covariance)[3])

        environmental_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_covariance))) # Pvalue


        # labile environmental coordination

        labile_environmental_covariance <- (totalResidualCV - pureResidualCV) / (sqrt( ((totalPhyloVar1 + totalNonPhyloVar1) * (totalPhyloVar2 + totalNonPhyloVar2)) ))

        # Significance
        labile_environmental_covariance_lCI <- as.numeric(RChronoModel::CredibleInterval(labile_environmental_covariance)[2])
        labile_environmental_covariance_hCI <- as.numeric(RChronoModel::CredibleInterval(labile_environmental_covariance)[3])

        labile_environmental_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(labile_environmental_covariance))) # Pvalue

        # total environmental coordination

        totalEnvironmentalCoordination <- environmental_phylogenetic_covariance + labile_environmental_covariance

        # Significance
        totalEnvironmentalCoordination_lCI <- as.numeric(RChronoModel::CredibleInterval(totalEnvironmentalCoordination)[2])
        totalEnvironmentalCoordination_hCI <- as.numeric(RChronoModel::CredibleInterval(totalEnvironmentalCoordination)[3])

        totalEnvironmentalCoordination_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(totalEnvironmentalCoordination))) # Pvalue

        # residual coordination

        residual_covariance <- total_covariance - (non_attributed_phylogenetic_covariance + environmental_phylogenetic_covariance + labile_environmental_covariance)

        # Significance
        residual_covariance_lCI <- as.numeric(RChronoModel::CredibleInterval(residual_covariance)[2])
        residual_covariance_hCI <- as.numeric(RChronoModel::CredibleInterval(residual_covariance)[3])

        residual_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(residual_covariance))) # Pvalue

        # results

        CVPartitionResults$covariancePartition <-  cbind(CVPartitionResults$covariancePartition,
                                                                 "environmental_variable" = paste0(environmental_variables, collapse = ", "),
                                                                 "non_attributed_phylogenetic_correlation" = mean(non_attributed_phylogenetic_covariance),
                                                                 "environmental_phylogenetic_correlation" = mean(environmental_phylogenetic_covariance),
                                                                 # "Total_environmental_coordination" = mean(totalEnvironmentalCoordination),
                                                                 "labile_environmental_correlation" = mean(labile_environmental_covariance),
                                                                 "residual_correlation" = mean(residual_covariance),

                                                                 "p_value_non_attributed_phylogenetic_correlation" = non_attributed_phylogenetic_covariance_pvalue,
                                                                 "low_CI_non_attributed_phylogenetic_correlation" = non_attributed_phylogenetic_covariance_lCI,
                                                                 "high_CI_non_attributed_phylogenetic_correlation" = non_attributed_phylogenetic_covariance_hCI,
                                                                 "CI_significance_non_attributed_phylogenetic_correlation" = ifelse(non_attributed_phylogenetic_covariance_lCI < 0 && non_attributed_phylogenetic_covariance_hCI > 0, "no", "yes"),

                                                                 "p_value_environmental_phylogenetic_correlation" = environmental_phylogenetic_covariance_pvalue,
                                                                 "low_CI_environmental_phylogenetic_correlation" = environmental_phylogenetic_covariance_lCI,
                                                                 "high_CI_environmental_phylogenetic_correlation" = environmental_phylogenetic_covariance_hCI,
                                                                 "CI_significance_environmental_phylogenetic_correlation" = ifelse(environmental_phylogenetic_covariance_lCI < 0 && environmental_phylogenetic_covariance_hCI > 0, "no", "yes"),

                                                                 # "Pvalue_Total_environmental_coordination" = totalEnvironmentalCoordination_pvalue,

                                                                 "p_value_labile_environmental_correlation" = labile_environmental_covariance_pvalue,
                                                                 "low_CI_labile_environmental_correlation" = labile_environmental_covariance_lCI,
                                                                 "high_CI_labile_environmental_correlation" = labile_environmental_covariance_hCI,
                                                                 "CI_significance_labile_environmental_correlation" = ifelse(labile_environmental_covariance_lCI < 0 && labile_environmental_covariance_hCI > 0, "no", "yes"),

                                                                 "p_value_residual_correlation" = residual_covariance_pvalue,
                                                                 "low_CI_residual_correlation" = residual_covariance_lCI,
                                                                 "high_CI_residual_correlation" = residual_covariance_hCI,
                                                                 "CI_significance_residual_correlation" = ifelse(residual_covariance_lCI < 0 && residual_covariance_hCI > 0, "no", "yes")
        )

        CVPartitionResults$covariancePartitionDistributions[["non_attributed_phylogenetic_correlation"]] <- non_attributed_phylogenetic_covariance
        CVPartitionResults$covariancePartitionDistributions[["environmental_phylogenetic_correlation"]] <- environmental_phylogenetic_covariance
        # CVPartitionResults$covariancePartitionDistributions[["totalEnvironmentalCoordination"]] <- totalEnvironmentalCoordination
        CVPartitionResults$covariancePartitionDistributions[["labile_environmental_correlation"]] <- labile_environmental_covariance
        CVPartitionResults$covariancePartitionDistributions[["residual_correlation"]] <- residual_covariance

        CVPartitionResults$modelPhyloEnv <- mdlPhyloEnv
        CVPartitionResults$model.diagnostics <- model.diagnostics

        # add to all traits results
        traitsCVPartitionResults$covariancePartition <- rbind(traitsCVPartitionResults$covariancePartition,
                                                                    CVPartitionResults$covariancePartition)

        traitsCVPartitionResults$model.diagnostics <- rbind(traitsCVPartitionResults$model.diagnostics,
                                                                     CVPartitionResults$model.diagnostics)
        traitsCVPartitionResults$individual.models.results[[model]] <- CVPartitionResults

      } # end evaluation if model already exists in results

    } # end bucle for all traits
  }

  # message("Model structure used:")
  # message(multi_mdls.str)
  # message("Phylogenetic signal results:")
  # message(CVPartitionResults$covariancePartition)

  # save results

  if(save){
    if(!is.null(environmental_variables)){
      assign(paste0("traitsCovariancePartitionResults_", environmental_variables, collapse = "_"), CVPartitionResults)
      save(list = paste0("traitsCovariancePartitionResults_", environmental_variables, collapse = "_"), file = results.file)
      message(results.file)
    } else{
      save(list = paste0("traitsCovariancePartitionResults"), file = results.file)
      message(results.file)
    }
  }
  return(traitsCVPartitionResults)
}
