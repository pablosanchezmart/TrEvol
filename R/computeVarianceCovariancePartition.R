#' Phylogenetic variance and covariance partition including (or not) environment
#'
#'Compute total variance and covariance and separate it into its phylogenetic and non-phylogenetic components.
#' If environmental variable is included, four variance and covariance is separated into four components: non-attributed phylogenetic variance-covariance, environmental
#'phylogenetic variance-covariance, non-phylogenetic (labile) environmental variance-covariance and residual variance-covariance.
#'
#'The function can automaticallys save results. For that, first run initializeTrEvol to create the folders and subfolders structure.
#'
#' @param traits (*character*). Name of the trait or list of traits to calculate variances and covariances. They must be in the dataset.
#' @param environmental_variable (*character*). Names of the environmental variables They must be contained in dataset.
#' @param dataset (*data frame*). Data frame containing the trait of interest and a column describing terminal appearing as tip labels of the phylogeny.
#' @param terminal_taxa (*character*). Terminal taxon as named in the dataset (e.g., species).
#' @param phylogeny (*phylo*) Phylogeny with tip labels contained in dataset$animal
#' @param model_specifications (*list*). Mcmcglmm models specifications as specified by the defineModelsSpecification of this package. If not defined, the function internally uses the default in defineModelsSpecification function.
#' @param show_relative_variance (*logical*). If true, variance is shown relative to the total variance (recommended when using different variables). Otherwise, absolute variances are reported.
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
#' # Compute variance-covariance structure for simulated traits using default parameters
#' variance_covariance_results <- computeVarianceCovariancePartition(
#' traits = c("phylo_G1_trait1", "phylo_G1_trait2"),
#' environmental_variable = "phylo_G1_env",
#' data = simulated_traits.data$data
#' phylogeny = simulated_traits.data$phylogeny
#' )
#' }
computeVarianceCovariancePartition <- function(traits = NULL,
                                               environmental_variable = NULL,
                                               dataset = NULL,
                                               terminal_taxa = NULL,
                                               phylogeny = NULL,
                                               model_specifications = NULL,
                                               show_relative_variance = T,
                                               force_run = T,
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

  # Results structure
  traitsVCVPartitionResults <- list()
  traitsVCVPartitionResults$varianceResults <- data.frame()
  traitsVCVPartitionResults$covarianceResults <- data.frame()
  traitsVCVPartitionResults$models.diagnostics <- data.frame()
  traitsVCVPartitionResults$individual.models.results <- list()


  ### MODELS STRUCTURE --------------------------------------------------------- ####

  multi_mdls.str <- expand.grid(traits, traits) # all possible pairwise combinations between variable
  names(multi_mdls.str) <- c("trait1", "trait2")
  multi_mdls.str$environment <- environmental_variable

  # This is needed in order to avoid models with the same response variables but in different order (which are equivalent models)
  multi_mdls.str$traits <- NA
  for(i in 1:length(multi_mdls.str$trait1)){
    both_traits <- c(as.character(multi_mdls.str[i, "trait1"]), as.character(multi_mdls.str[i, "trait2"]))
    both_traits <- sort(both_traits)
    multi_mdls.str$traits[i] <- paste0(both_traits[1], ", ", both_traits[2])
    multi_mdls.str$trait1[i] <- as.character(both_traits[[1]])
    multi_mdls.str$trait2[i] <-  as.character(both_traits[[2]])
  }

  if(!is.null(environmental_variable)){
    multi_mdls.str$type <- paste0("tri_", multi_mdls.str$trait1, "_", multi_mdls.str$trait2, "_", multi_mdls.str$environment)
    multi_mdls.str$fix.frml <- paste0("cbind(", multi_mdls.str$traits, ", ", multi_mdls.str$environment, ") ~ trait-1")
  } else{
    multi_mdls.str$type <- paste0("bi_", multi_mdls.str$trait1, "_", multi_mdls.str$trait2)
    multi_mdls.str$fix.frml <- paste0("cbind(", multi_mdls.str$traits, ") ~ trait-1")
  }
  multi_mdls.str$ran.frml <- "~ us(trait):animal"

  multi_mdls.str <- multi_mdls.str %>%
    dplyr::filter(!trait1 == trait2) %>%
    dplyr::filter(!duplicated(traits)) %>%
    dplyr::select(type, traits, trait1, trait2, fix.frml, ran.frml)

  ### RUN MODELS --------------------------------------------------------------- ####

  # if (file.exists(results.file) && isFALSE(force_run)) {
  #   print("loanding previous results")
  #   load(file = results.file)
  # }


  # run models and extract results
  for (model in multi_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(traitsVCVPartitionResults$individual.models.results) | force_run) {

      message(paste0("Running variance-covariance calculation: ", model))
      model.descr <- multi_mdls.str %>%
        dplyr::filter(type == model)

      trait1 <- as.character(model.descr$trait1)
      trait2 <- as.character(model.descr$trait2)

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait1, trait2, environmental_variable))

      if (is.null(model_specifications)) {
        warning("Using default model specificatios. Use defineModelsSpecifications() output on model_specifications argument to set them manually.")
        model_specifications <- defineModelsSpecifications()
      }

      # model
      if(!is.null(environmental_variable)){
        fix.frml <- paste0("cbind(", trait1, ", ", trait2, ", ", environmental_variable, ") ~ trait-1")

        mdlPhyloEnv <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                          random = ~ us(trait):animal, rcov = ~us(trait):units,
                                          data= modellingData$dta, pedigree = modellingData$phylo,
                                          family = c("gaussian", "gaussian", "gaussian"),
                                          prior = model_specifications$triresponse_prior,
                                          nitt = model_specifications$number_interations,
                                          burnin = model_specifications$burning_iterations,
                                          thin = model_specifications$thinning_iterations,
                                          verbose = F)

      } else{
        fix.frml <- paste0("cbind(", trait1, ", ", trait2, ") ~ trait-1")

        mdlPhyloEnv <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                          random = ~ us(trait):animal, rcov = ~us(trait):units,
                                          data= modellingData$dta, pedigree = modellingData$phylo,
                                          family = c("gaussian", "gaussian"),
                                          prior = model_specifications$biresponse_prior,
                                          nitt = model_specifications$number_interations,
                                          burnin = model_specifications$burning_iterations,
                                          thin = model_specifications$thinning_iterations,
                                          verbose = F)

      }

      mdlPhyloEnv$name <- fix.frml

      model.diagnostics <- diagnoseModels(model = mdlPhyloEnv, verbose = verbose)

      #### VARIANCE CALCULATIONS ----------------------------------------------- ####

      ### Total variance ####

      # Trait 1

      total_variance_t1 <- abs(mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".animal")] +
                                 mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")])

      total_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_variance_t1))) # Pvalue

      # Trait 2

      total_variance_t2 <- abs(mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".animal")]) +
        abs(mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".units")])

      total_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_variance_t2))) # Pvalue


      ### Total phylogenetic variance ####

      # Trait 1

      total_phylogenetic_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".animal")]

      if(show_relative_variance){
        total_phylogenetic_variance_t1 <- total_phylogenetic_variance_t1 / total_variance_t1
      }

      total_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_variance_t1))) # Pvalue

      # Trait 2

      total_phylogenetic_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".animal")]

      if(show_relative_variance){
        total_phylogenetic_variance_t2 <- total_phylogenetic_variance_t2 / total_variance_t2
      }

      total_phylogenetic_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_variance_t2))) # Pvalue


      ### Total non-phylogenetic variance ####

      # Trait 1

      total_non_phylogenetic_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")]

      if(show_relative_variance){
        total_non_phylogenetic_variance_t1 <- total_non_phylogenetic_variance_t1 / total_variance_t1
      }

      total_non_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_variance_t1))) # Pvalue

      # Trait 2

      total_non_phylogenetic_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait1, ".units")]

      if(show_relative_variance){
        total_non_phylogenetic_variance_t2 <- total_non_phylogenetic_variance_t2 / total_variance_t2
      }

      total_non_phylogenetic_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_variance_t2))) # Pvalue


      if(!is.null(environmental_variable)){


        ### Non-attributed phylogenetic variance ####

        # Trait 1

        non_attributed_phylogenetic_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".animal")] -     # VAR(u1)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".animal")] /                  # COV(UE,u1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]  # VAR(uE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]      # VAR(uE)

        if(show_relative_variance){
          non_attributed_phylogenetic_variance_t1 <- non_attributed_phylogenetic_variance_t1 / total_variance_t1
        }

        non_attributed_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(non_attributed_phylogenetic_variance_t1))) # Pvalue

        # Trait 2

        non_attributed_phylogenetic_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".animal")] -     # VAR(u1)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".animal")] /                  # COV(UE,u1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]  # VAR(uE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]      # VAR(uE)

        if(show_relative_variance){
          non_attributed_phylogenetic_variance_t2 <- non_attributed_phylogenetic_variance_t2 / total_variance_t2
        }

        non_attributed_phylogenetic_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(non_attributed_phylogenetic_variance_t2))) # Pvalue


        ### Environmental phylogenetic variance ####

        # Trait 1

        environmental_phylogenetic_variance_t1 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".animal")] /                # COV(UE,u1)
            mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]  # VAR(uE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]  # VAR(uE)

        if(show_relative_variance){
          environmental_phylogenetic_variance_t1 <- environmental_phylogenetic_variance_t1 / total_variance_t1
        }

        environmental_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_variance_t1))) # Pvalue

        # Trait 2

        environmental_phylogenetic_variance_t2 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".animal")] /                  # COV(UE,u2)
            mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]  # VAR(uE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")]    # VAR(uE)

        if(show_relative_variance){
          environmental_phylogenetic_variance_t2 <- environmental_phylogenetic_variance_t2 / total_variance_t2
        }

        environmental_phylogenetic_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_variance_t2))) # Pvalue


        ### Labile (non-phylogenetic) environmental variance ####

        # Trait 1

        labile_environmental_variance_t1 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".units")] /                  # COV(eE,e1)
            mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]  # VAR(eE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]    # VAR(eE)

        if(show_relative_variance){
          labile_environmental_variance_t1 <- labile_environmental_variance_t1 / total_variance_t1
        }

        labile_environmental_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(labile_environmental_variance_t1))) # Pvalue

        # Trait 2

        labile_environmental_variance_t2 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".units")] /                  # COV(eE,e2)
            mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]  # VAR(eE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]    # VAR(eE)

        if(show_relative_variance){
          labile_environmental_variance_t2 <- labile_environmental_variance_t2 / total_variance_t2
        }

        labile_environmental_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(labile_environmental_variance_t2))) # Pvalue


        ### Residual variance ####

        # Trait 1

        residual_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")] -         # VAR(e1)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".units")] /                  # COV(eE,e1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]  # VAR(eE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]      # VAR(eE)

        if(show_relative_variance){
          residual_variance_t1 <- residual_variance_t1 / total_variance_t1
        }

        residual_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(residual_variance_t1))) # Pvalue

        # Trait 2

        residual_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".units")] -         # VAR(e2)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".units")] /                  # COV(eE,e1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]  # VAR(eE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")]      # VAR(eE)

        if(show_relative_variance){
          residual_variance_t2 <- residual_variance_t2 / total_variance_t2
        }

        residual_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(residual_variance_t2))) # Pvalue
      }


      if(length(traits) >= 2){
        #### COVARIANCE CALCULATIONS --------------------------------------------- ####
        ### Total covariance ####

        total_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")] +                 # COV(u1, u2)
          mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]                                      # COV(e1, e2)

        # Correlation
          total_covariance <- total_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )

        total_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_covariance))) # Pvalue


        ### Total phylogenetic covariance ####

        total_phylogenetic_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")]                 # COV(u1, u2)

        # Correlation
          total_phylogenetic_covariance <- total_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )

        total_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_covariance))) # Pvalue


        ### Total non-phylogenetic covariance ####

        total_non_phylogenetic_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]                 # COV(u1, u2)

        # Correlation
          total_non_phylogenetic_covariance <- total_non_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )

        total_non_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_covariance))) # Pvalue


        if(!is.null(environmental_variable)){


          ### Non-attributed phylogenetic covariance ####

          non_attributed_phylogenetic_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")] -                               # COV(u1, u2)
            (
              (
                mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".animal")] *               # COV(u1, uE)
                  mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".animal")]                 # COV(u2, uE)
              )/
                mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")] # var(uE)
            )

          # Correlation
            non_attributed_phylogenetic_covariance <- non_attributed_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )

          non_attributed_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(non_attributed_phylogenetic_covariance))) # Pvalue


          ### Environmental phylogenetic covariance ####

          environmental_phylogenetic_covariance <- (
            mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".animal")] *               # COV(u1, uE)
              mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".animal")]                 # COV(u2, uE)
          )/
            mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".animal")] # var(uE)

          # Correlation
            environmental_phylogenetic_covariance <- environmental_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )

          environmental_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_covariance))) # Pvalue


          ### Labile (non-phylogenetic) environmental covariance ####

          labile_environmental_covariance <- (
            mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".units")] *               # COV(e1, eE)
              mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".units")]                 # COV(e2, eE)
          )/
            mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")] # var(eE)

          # Correlation
            labile_environmental_covariance <- labile_environmental_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )

          labile_environmental_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(labile_environmental_covariance))) # Pvalue


          ### Residual covariance ####

          residual_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")] -                               # COV(e1, e2)
            (
              (
                mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental_variable, ".units")] *               # COV(e1, eE)
                  mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental_variable, ".units")]                 # COV(e2, eE)
              )/
                mdlPhyloEnv$VCV[, paste0("trait", environmental_variable, ":trait", environmental_variable, ".units")] # var(eE)
            )

          # Correlation
            residual_covariance <- residual_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )

          residual_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(residual_covariance))) # Pvalue
        }



      }
      ### RESULTS -------------------------------------------------------------- ####

      ### Variance results ####

      VCVPartitionResults <- list()

      if(!trait1 %in% traitsVCVPartitionResults$varianceResults$Trait){

        # trait 1

        VCVPartitionResults$variancePartition <- data.frame("trait" = trait1,
                                                            "number_observations" = length(modellingData$dta$animal),
                                                            "phylogenetic_variance" = mean(total_phylogenetic_variance_t1),
                                                            "non_phylogenetic_variance" = mean(total_non_phylogenetic_variance_t1),
                                                            "p_value_phylogenetic_variance" = total_phylogenetic_variance_t1_pvalue
        )

        VCVPartitionResults$variancePartitionDistributions[[paste0("phylogenetic_variance_", trait1)]] <- total_phylogenetic_variance_t1
        VCVPartitionResults$variancePartitionDistributions[[paste0("non_phylogenetic_variance_", trait1)]] <- total_non_phylogenetic_variance_t1

        if(!is.null(environmental_variable)){
          VCVPartitionResults$variancePartition <- cbind(VCVPartitionResults$variancePartition,
                                                         "environmental_variable" = paste0(environmental_variable, collapse = ", "),
                                                         "non_attributed_phylogenetic_variance" = mean(non_attributed_phylogenetic_variance_t1),
                                                         "environmental_phylogenetic_variance" = mean(environmental_phylogenetic_variance_t1),
                                                         "labile_environmental_variance" = mean(labile_environmental_variance_t1),
                                                         "residual_variance" = mean(residual_variance_t1),
                                                         "p_value_non_attributed_phylogenetic_variance" = non_attributed_phylogenetic_variance_t1_pvalue,
                                                         "p_value_environmental_phylogenetic_variance" = environmental_phylogenetic_variance_t1_pvalue,
                                                         "p_value_labile_environmental_variance" = environmental_phylogenetic_variance_t1_pvalue
          )

          VCVPartitionResults$variancePartitionDistributions[[paste0("non_attributed_phylogenetic_variance", trait1)]] <- non_attributed_phylogenetic_variance_t1
          VCVPartitionResults$variancePartitionDistributions[[paste0("environmental_phylogenetic_variance", trait1)]] <- environmental_phylogenetic_variance_t1
          VCVPartitionResults$variancePartitionDistributions[[paste0("labile_environmental_variance", trait1)]] <- labile_environmental_variance_t1
          VCVPartitionResults$variancePartitionDistributions[[paste0("residual_variance", trait1)]] <- residual_variance_t1
        }
      }

      if(!trait2 %in% traitsVCVPartitionResults$varianceResults$trait){

        # trait 2

        variancePartition_t2 <- data.frame("trait" = trait2,
                                           "number_observations" = length(modellingData$dta$animal),
                                           "phylogenetic_variance" = mean(total_phylogenetic_variance_t2),
                                           "non_phylogenetic_variance" = mean(total_non_phylogenetic_variance_t2),
                                           "p_value_phylogenetic_variance" = total_phylogenetic_variance_t2_pvalue
        )

        VCVPartitionResults$variancePartitionDistributions[[paste0("phylogenetic_variance_", trait2)]] <- total_phylogenetic_variance_t2
        VCVPartitionResults$variancePartitionDistributions[[paste0("non_phylogenetic_variance_", trait2)]] <- total_non_phylogenetic_variance_t2

        if(!is.null(environmental_variable)){
          variancePartition_t2 <- cbind(variancePartition_t2,
                                        "environmental_variable" = paste0(environmental_variable, collapse = ", "),
                                        "non_attributed_phylogenetic_variance" = mean(non_attributed_phylogenetic_variance_t2),
                                        "environmental_phylogenetic_variance" = mean(environmental_phylogenetic_variance_t2),
                                        "labile_environmental_variance" = mean(labile_environmental_variance_t2),
                                        "residual_variance" = mean(residual_variance_t2),
                                        "p_value_non_attributed_phylogenetic_variance" = non_attributed_phylogenetic_variance_t2_pvalue,
                                        "p_value_environmental_phylogenetic_variance" = environmental_phylogenetic_variance_t2_pvalue,
                                        "p_value_labile_environmental_variance" = environmental_phylogenetic_variance_t2_pvalue
          )

          VCVPartitionResults$variancePartitionDistributions[[paste0("non_attributed_phylogenetic_variance", trait2)]] <- non_attributed_phylogenetic_variance_t2
          VCVPartitionResults$variancePartitionDistributions[[paste0("environmental_phylogenetic_variance", trait2)]] <- environmental_phylogenetic_variance_t2
          VCVPartitionResults$variancePartitionDistributions[[paste0("labile_environmental_variance", trait2)]] <- labile_environmental_variance_t2
          VCVPartitionResults$variancePartitionDistributions[[paste0("residual_variance", trait2)]] <- residual_variance_t2
        }

        if(!trait1 %in% traitsVCVPartitionResults$varianceResults$trait){
          VCVPartitionResults$variancePartition <- rbind(VCVPartitionResults$variancePartition, variancePartition_t2)
        } else {
          VCVPartitionResults$variancePartition <- variancePartition_t2
        }

      }


      if(length(traits) >= 2){

        ### Covariance results ####

        VCVPartitionResults$covariancePartition <-  data.frame("trait_1" = trait1,
                                                               "trait_2" = trait2,
                                                               "number_observations" = length(modellingData$dta$animal),

                                                               "total_correlation" = mean(total_covariance),
                                                               "phylogenetic_correlation" = mean(total_phylogenetic_covariance),
                                                               "non_phylogenetic_correlation" = mean(total_non_phylogenetic_covariance),
                                                               "p_value_total_correlation" = total_covariance_pvalue,
                                                               "p_value_phylogenetic_correlation" = total_phylogenetic_covariance_pvalue,
                                                               "p_value_non_phylogenetic_correlation" = total_non_phylogenetic_covariance_pvalue
        )

        VCVPartitionResults$covariancePartitionDistributions <- list("total_correlation" = total_covariance,
                                                                     "phylogenetic_correlation" = total_phylogenetic_covariance,
                                                                     "non_phylogenetic_correlation" = total_non_phylogenetic_covariance)


        if(!is.null(environmental_variable)){
          VCVPartitionResults$covariancePartition <-  cbind(VCVPartitionResults$covariancePartition,
                                                            "environmental_variable" = paste0(environmental_variable, collapse = ", "),
                                                            "non_attributed_phylogenetic_correlation" = mean(non_attributed_phylogenetic_covariance),
                                                            "environmental_phylogenetic_correlation" = mean(environmental_phylogenetic_covariance),
                                                            "labile_environmental_correlation" = mean(labile_environmental_covariance),
                                                            "residual_correlation" = mean(residual_covariance),

                                                            "p_value_non_attributed_phylogenetic_correlation" = non_attributed_phylogenetic_covariance_pvalue,
                                                            "p_value_environmental_phylogenetic_correlation" = environmental_phylogenetic_covariance_pvalue,
                                                            "p_value_labile_environmental_correlation" = labile_environmental_covariance_pvalue,
                                                            "p_value_residual_correlation" = residual_covariance_pvalue
          )

          VCVPartitionResults$covariancePartitionDistributions[["non_attributed_phylogenetic_correlation"]] <- non_attributed_phylogenetic_covariance
          VCVPartitionResults$covariancePartitionDistributions[["environmental_phylogenetic_correlation"]] <- environmental_phylogenetic_covariance
          VCVPartitionResults$covariancePartitionDistributions[["labile_environmental_correlation"]] <- labile_environmental_covariance
          VCVPartitionResults$covariancePartitionDistributions[["residual_correlation"]] <- residual_covariance
        }
      }

        VCVPartitionResults$modelPhyloEnv <- mdlPhyloEnv
        VCVPartitionResults$model.diagnostics <- model.diagnostics

        # add to all traits results
        traitsVCVPartitionResults$covarianceResults <- rbind(traitsVCVPartitionResults$covarianceResults,
                                                             VCVPartitionResults$covariancePartition)

        traitsVCVPartitionResults$varianceResults <-  rbind(traitsVCVPartitionResults$varianceResults,
                                                            VCVPartitionResults$variancePartition)

        # Calculate mean per trait in variance results

        traitsVCVPartitionResults$varianceResults <- stats::aggregate(traitsVCVPartitionResults$varianceResults[, -1],
                                                               by = list(traitsVCVPartitionResults$varianceResults$trait),
                                                               FUN = meanOrMode) %>%
          dplyr::rename(trait = Group.1)

        if(!is.null(environmental_variable)){
          traitsVCVPartitionResults$varianceResults$environmental_variable <- environmental_variable
        }

        traitsVCVPartitionResults$models.diagnostics <- rbind(traitsVCVPartitionResults$models.diagnostics,
                                                              VCVPartitionResults$model.diagnostics)
        traitsVCVPartitionResults$individual.models.results[[model]] <- VCVPartitionResults


    } # end evaluation if model already exists in results

  } # end bucle for all traits

  # message("Model structure used:")
  # message(multi_mdls.str)
  # print("Variance results:")
  # print(traitsVCVPartitionResults$varianceResults)
  # print("Covariance results:")
  # print(traitsVCVPartitionResults$covarianceResults)

  # Variance

  if(save){

    # save results
    results.file <- paste0(outputs.dir, "/models_outputs/")

    if(!is.null(environmental_variable)){
      assign(paste0("traits_variance_partition_results_", environmental_variable), traitsVCVPartitionResults$varianceResults)

      save(list = paste0("traits_variance_partition_results_", environmental_variable),
           file = paste0(results.file, paste0("traits_variance_partition_results_", environmental_variable, ".RData")))
      message(results.file)
    } else{
      assign(paste0("traits_variance_partition_results"), traitsVCVPartitionResults$varianceResults)

      save(list = paste0("traits_variance_partition_results"),
           file = paste0(results.file, "traits_variance_partition_results.RData"))
      message(results.file)
    }

  # Covariance

    if(!is.null(environmental_variable)){
      assign(paste0("traits_covariance_partition_results_", environmental_variable), traitsVCVPartitionResults$covarianceResults)

      save(list = paste0("traits_covariance_partition_results_", environmental_variable),
           file = paste0(results.file, paste0("traits_covariance_partition_results_", environmental_variable, ".RData")))
      message(results.file)
    } else{
      assign(paste0("traits_covariance_partition_results"), traitsVCVPartitionResults$covarianceResults)

      save(list = paste0("traits_covariance_partition_results"),
           file = paste0(results.file, "traits_covariance_partition_results.RData"))
      message(results.file)
    }
  }

  return(traitsVCVPartitionResults)
}
