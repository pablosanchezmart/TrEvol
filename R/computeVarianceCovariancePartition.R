
#' Variance partition including (or not) environment
#'
#' @param traits (character) Name of the trait or list of traits. It  must be contained in the dataset.
#' @param environmental.variable (character) Names of the environmental variables They must be contained in dataset.
#' @param dataset (data frame) Dataset containing the trait of interest and a column named "animal" describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param model.specifications (list) Mcmcglmm models specifications. See defineModelsSpecification.
#' @showCorrelations (logical) If true, relative variances and covariances are reported.
#' @param force.run (logical) If false, models already run are not runned again.
#' @param save (logical) If false, resulta re not saved.
#'
#' @return
#' @export
#'
#' @examples
computeVarianceCovariancePartition <- function(traits, environmental.variable = NULL, dataset, phylogeny, model.specifications = NULL, showCorrelations = T,
                                               force.run = T, save = T) {

  traitsVCVPartitionResults <- list()
  traitsVCVPartitionResults$varianceResults <- data.frame()
  traitsVCVPartitionResults$covarianceResults <- data.frame()
  traitsVCVPartitionResults$models.diagnostics <- data.frame()
  traitsVCVPartitionResults$individual.models.results <- list()

  ### MODELS STRUCTURE ####

  multi_mdls.str <- expand.grid(traits, traits) # all possible pairwise combinations between variable
  names(multi_mdls.str) <- c("trait1", "trait2")
  multi_mdls.str$environment <- environmental.variable

  # This is needed in order to avoid models with the same response variables but in different order (which are equivalent models)
  multi_mdls.str$traits <- NA
  for(i in 1:length(multi_mdls.str$trait1)){
    both_traits <- c(as.character(multi_mdls.str[i, "trait1"]), as.character(multi_mdls.str[i, "trait2"]))
    both_traits <- sort(both_traits)
    multi_mdls.str$traits[i] <- paste0(both_traits[1], ", ", both_traits[2])
    multi_mdls.str$trait1[i] <- as.character(both_traits[[1]])
    multi_mdls.str$trait2[i] <-  as.character(both_traits[[2]])
  }

  if(!is.null(environmental.variable)){
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

  ### RUN MODELS ####

  if(!is.null(environmental.variable)){
    results.file <- paste0(outputs.dir, "/models_outputs/traitsVCVPartitionResults", environmental.variable, ".RData")
  } else {
    results.file <- paste0(outputs.dir, "/models_outputs/traitsVCVPartitionResults.RData")
  }

  if (file.exists(results.file) && isFALSE(force.run)) {
    print("loanding previous results")
    load(file = results.file)
  }


  # run models and extract results
  for (model in multi_mdls.str$type) {

    # avoid running models already present in results
    if (!model %in% names(traitsVCVPartitionResults$individual.models.results) | force.run) {

      print(paste0("Running covariance calculation: ", model))
      model.descr <- multi_mdls.str %>%
        dplyr::filter(type == model)

      trait1 <- as.character(model.descr$trait1)
      trait2 <- as.character(model.descr$trait2)

      modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(trait1, trait2, environmental.variable))

      if (is.null(model.specifications)) {
        print("Using default model specificatios. Use defineModelsSpecifications() output on model.specifications argument to set them manually.")
        model.specifications <- defineModelsSpecifications()
      }

      # model
      if(!is.null(environmental.variable)){
        fix.frml <- paste0("cbind(", trait1, ", ", trait2, ", ", environmental.variable, ") ~ trait-1")

        mdlPhyloEnv <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                          random = ~ us(trait):animal, rcov = ~us(trait):units,
                                          data= modellingData$dta, pedigree = modellingData$phylo,
                                          family = c("gaussian", "gaussian", "gaussian"),
                                          prior = model.specifications$triresponse_prior,
                                          nitt = model.specifications$number_interations,
                                          burnin = model.specifications$burning_iterations,
                                          thin = model.specifications$thinning_iterations,
                                          verbose = F)

      } else{
        fix.frml <- paste0("cbind(", trait1, ", ", trait2, ") ~ trait-1")

        mdlPhyloEnv <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                                          random = ~ us(trait):animal, rcov = ~us(trait):units,
                                          data= modellingData$dta, pedigree = modellingData$phylo,
                                          family = c("gaussian", "gaussian"),
                                          prior = model.specifications$multiresponse_prior,
                                          nitt = model.specifications$number_interations,
                                          burnin = model.specifications$burning_iterations,
                                          thin = model.specifications$thinning_iterations,
                                          verbose = F)

      }

      mdlPhyloEnv$name <- fix.frml

      model.diagnostics <- diagnoseModels(model = mdlPhyloEnv)

      ### VARIANCE CALCULATIONS ####

      ## Total variance

      # Trait 1

      total_variance_t1 <- abs(mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".animal")] +
                                 mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")])

      total_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_variance_t1))) # Pvalue

      # Trait 2

      total_variance_t2 <- abs(mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".animal")]) +
        abs(mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".units")])

      total_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_variance_t2))) # Pvalue


      ## Total phylogenetic variance

      # Trait 1

      total_phylogenetic_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".animal")]

      if(showCorrelations){
        total_phylogenetic_variance_t1 <- total_phylogenetic_variance_t1 / total_variance_t1
      }

      total_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_variance_t1))) # Pvalue

      # Trait 2

      total_phylogenetic_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".animal")]

      if(showCorrelations){
        total_phylogenetic_variance_t2 <- total_phylogenetic_variance_t2 / total_variance_t2
      }

      total_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_variance_t2))) # Pvalue


      ## Total non-phylogenetic variance

      # Trait 1

      total_non_phylogenetic_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")]

      if(showCorrelations){
        total_non_phylogenetic_variance_t1 <- total_non_phylogenetic_variance_t1 / total_variance_t1
      }

      total_non_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_variance_t1))) # Pvalue

      # Trait 2

      total_non_phylogenetic_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait1, ".units")]

      if(showCorrelations){
        total_non_phylogenetic_variance_t2 <- total_non_phylogenetic_variance_t2 / total_variance_t2
      }

      total_non_phylogenetic_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_variance_t2))) # Pvalue


      if(!is.null(environmental.variable)){


        ## Pure phylogenetic variance

        # Trait 1

        pure_phylogenetic_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".animal")] -     # VAR(u1)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".animal")] /                  # COV(UE,u1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]  # VAR(uE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]      # VAR(uE)

        if(showCorrelations){
          pure_phylogenetic_variance_t1 <- pure_phylogenetic_variance_t1 / total_variance_t1
        }

        pure_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_phylogenetic_variance_t1))) # Pvalue

        # Trait 2

        pure_phylogenetic_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".animal")] -     # VAR(u1)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".animal")] /                  # COV(UE,u1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]  # VAR(uE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]      # VAR(uE)

        if(showCorrelations){
          pure_phylogenetic_variance_t2 <- pure_phylogenetic_variance_t2 / total_variance_t2
        }

        pure_phylogenetic_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_phylogenetic_variance_t2))) # Pvalue


        ## Environmental phylogenetic variance

        # Trait 1

        environmental_phylogenetic_variance_t1 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".animal")] /                # COV(UE,u1)
            mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]  # VAR(uE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]  # VAR(uE)

        if(showCorrelations){
          environmental_phylogenetic_variance_t1 <- environmental_phylogenetic_variance_t1 / total_variance_t1
        }

        environmental_phylogenetic_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_variance_t1))) # Pvalue

        # Trait 2

        environmental_phylogenetic_variance_t2 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".animal")] /                  # COV(UE,u2)
            mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]  # VAR(uE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")]    # VAR(uE)

        if(showCorrelations){
          environmental_phylogenetic_variance_t2 <- environmental_phylogenetic_variance_t2 / total_variance_t2
        }

        environmental_phylogenetic_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_variance_t2))) # Pvalue


        ## Pure environmental variance

        # Trait 1

        pure_environmental_variance_t1 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".units")] /                  # COV(eE,e1)
            mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]  # VAR(eE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]    # VAR(eE)

        if(showCorrelations){
          pure_environmental_variance_t1 <- pure_environmental_variance_t1 / total_variance_t1
        }

        pure_environmental_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_environmental_variance_t1))) # Pvalue

        # Trait 2

        pure_environmental_variance_t2 <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".units")] /                  # COV(eE,e2)
            mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]  # VAR(eE)
        )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]    # VAR(eE)

        if(showCorrelations){
          pure_environmental_variance_t2 <- pure_environmental_variance_t2 / total_variance_t2
        }

        pure_environmental_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_environmental_variance_t2))) # Pvalue


        ## Pure residual variance

        # Trait 1

        pure_residual_variance_t1 <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait1, ".units")] -         # VAR(e1)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".units")] /                  # COV(eE,e1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]  # VAR(eE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]      # VAR(eE)

        if(showCorrelations){
          pure_residual_variance_t1 <- pure_residual_variance_t1 / total_variance_t1
        }

        pure_residual_variance_t1_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_residual_variance_t1))) # Pvalue

        # Trait 2

        pure_residual_variance_t2 <- mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", trait2, ".units")] -         # VAR(e2)
          (
            mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".units")] /                  # COV(eE,e1)
              mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]  # VAR(eE)
          )^2 *
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")]      # VAR(eE)

        if(showCorrelations){
          pure_residual_variance_t2 <- pure_residual_variance_t2 / total_variance_t2
        }

        pure_residual_variance_t2_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_residual_variance_t2))) # Pvalue
      }


      ### COVARIANCE CALCULATIONS ####

      ## Total covariance

      total_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")] +                 # COV(u1, u2)
        mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]                                      # COV(e1, e2)

      if(showCorrelations){
        total_covariance <- total_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )
      }

      total_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_covariance))) # Pvalue

      ## Total phylogenetic covariance

      total_phylogenetic_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")]                 # COV(u1, u2)

      if(showCorrelations){
        total_phylogenetic_covariance <- total_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )
      }

      total_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_phylogenetic_covariance))) # Pvalue

      ## Total non-phylogenetic covariance

      total_non_phylogenetic_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")]                 # COV(u1, u2)

      if(showCorrelations){
        total_non_phylogenetic_covariance <- total_non_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )
      }

      total_non_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(total_non_phylogenetic_covariance))) # Pvalue


      if(!is.null(environmental.variable)){

        ## Pure phylogenetic covariance

        pure_phylogenetic_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".animal")] -                               # COV(u1, u2)
          (
            (
              mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".animal")] *               # COV(u1, uE)
                mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".animal")]                 # COV(u2, uE)
            )/
              mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")] # var(uE)
          )

        if(showCorrelations){
          pure_phylogenetic_covariance <- pure_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )
        }

        pure_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_phylogenetic_covariance))) # Pvalue


        ## Environmental phylogenetic covariance

        environmental_phylogenetic_covariance <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".animal")] *               # COV(u1, uE)
            mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".animal")]                 # COV(u2, uE)
        )/
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".animal")] # var(uE)

        if(showCorrelations){
          environmental_phylogenetic_covariance <- environmental_phylogenetic_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )
        }

        environmental_phylogenetic_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(environmental_phylogenetic_covariance))) # Pvalue


        ## Pure environmental covariance

        pure_environmental_covariance <- (
          mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".units")] *               # COV(e1, eE)
            mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".units")]                 # COV(e2, eE)
        )/
          mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")] # var(eE)

        if(showCorrelations){
          pure_environmental_covariance <- pure_environmental_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )
        }

        pure_environmental_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_environmental_covariance))) # Pvalue


        ## Pure residual covariance

        pure_residual_covariance <- mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", trait2, ".units")] -                               # COV(e1, e2)
          (
            (
              mdlPhyloEnv$VCV[, paste0("trait", trait1, ":trait", environmental.variable, ".units")] *               # COV(e1, eE)
                mdlPhyloEnv$VCV[, paste0("trait", trait2, ":trait", environmental.variable, ".units")]                 # COV(e2, eE)
            )/
              mdlPhyloEnv$VCV[, paste0("trait", environmental.variable, ":trait", environmental.variable, ".units")] # var(eE)
          )

        if(showCorrelations){
          pure_residual_covariance <- pure_residual_covariance  / (sqrt( ((total_variance_t1) * (total_variance_t2)) ) )
        }

        pure_residual_covariance_pvalue <- 2*(1 - as.numeric(bayestestR::p_direction(pure_residual_covariance))) # Pvalue
      }


      ### RESULTS ####

      ## Variance results

      VCVPartitionResults <- list()

      if(!trait1 %in% traitsVCVPartitionResults$varianceResults$Trait){

        # Trait 1

        VCVPartitionResults$variancePartition <- data.frame("Trait" = trait1,
                                                            "N" = length(modellingData$dta$animal),
                                                            "Total_phylogenetic_conservatism" = mean(total_phylogenetic_variance_t1),
                                                            "Total_non_phylogenetic" = mean(total_non_phylogenetic_variance_t1)
        )

        VCVPartitionResults$variancePartitionDistributions[[paste0("Total_phylogenetic_conservatism_", trait1)]] <- total_phylogenetic_variance_t1
        VCVPartitionResults$variancePartitionDistributions[[paste0("Total_non_phylogenetic_", trait1)]] <- total_non_phylogenetic_variance_t1

        if(!is.null(environmental.variable)){
          VCVPartitionResults$variancePartition <- cbind(VCVPartitionResults$variancePartition,
                                                         "Environmental_variables" = paste0(environmental.variable, collapse = ", "),
                                                         "Pure_phylogenetic_conservatism" = mean(pure_phylogenetic_variance_t1),
                                                         "Phylogenetic_niche_conservatism" = mean(environmental_phylogenetic_variance_t1),
                                                         "Pure_environmental" = mean(pure_environmental_variance_t1),
                                                         "Residual" = mean(pure_residual_variance_t1)
          )

          VCVPartitionResults$variancePartitionDistributions[[paste0("Pure_phylogenetic_conservatism", trait1)]] <- pure_phylogenetic_variance_t1
          VCVPartitionResults$variancePartitionDistributions[[paste0("Phylogenetic_niche_conservatism", trait1)]] <- environmental_phylogenetic_variance_t1
          VCVPartitionResults$variancePartitionDistributions[[paste0("pure_environmental", trait1)]] <- pure_environmental_variance_t1
          VCVPartitionResults$variancePartitionDistributions[[paste0("residual", trait1)]] <- pure_residual_variance_t1
        }
      }

      if(!trait2 %in% traitsVCVPartitionResults$varianceResults$Trait){

        # Trait 2

        variancePartition_t2 <- data.frame("Trait" = trait2,
                                           "N" = length(modellingData$dta$animal),
                                           "Total_phylogenetic_conservatism" = mean(total_phylogenetic_variance_t2),
                                           "Total_non_phylogenetic" = mean(total_non_phylogenetic_variance_t2)
        )

        VCVPartitionResults$variancePartitionDistributions[[paste0("Total_phylogenetic_conservatism_", trait2)]] <- total_phylogenetic_variance_t2
        VCVPartitionResults$variancePartitionDistributions[[paste0("Total_non_phylogenetic_", trait2)]] <- total_non_phylogenetic_variance_t2

        if(!is.null(environmental.variable)){
          variancePartition_t2 <- cbind(variancePartition_t2,
                                        "Environmental_variables" = paste0(environmental.variable, collapse = ", "),
                                        "Pure_phylogenetic_conservatism" = mean(pure_phylogenetic_variance_t2),
                                        "Phylogenetic_niche_conservatism" = mean(environmental_phylogenetic_variance_t2),
                                        "Pure_environmental" = mean(pure_environmental_variance_t2),
                                        "Residual" = mean(pure_residual_variance_t2)
          )

          VCVPartitionResults$variancePartitionDistributions[[paste0("Pure_phylogenetic_conservatism_", trait2)]] <- pure_phylogenetic_variance_t2
          VCVPartitionResults$variancePartitionDistributions[[paste0("Phylogenetic_niche_conservatism", trait2)]] <- environmental_phylogenetic_variance_t2
          VCVPartitionResults$variancePartitionDistributions[[paste0("pure_environmental", trait2)]] <- pure_environmental_variance_t2
          VCVPartitionResults$variancePartitionDistributions[[paste0("residual", trait2)]] <- pure_residual_variance_t2
        }

        if(!trait1 %in% traitsVCVPartitionResults$varianceResults$Trait){
          VCVPartitionResults$variancePartition <- rbind(VCVPartitionResults$variancePartition, variancePartition_t2)
        } else {
          VCVPartitionResults$variancePartition <- variancePartition_t2
        }

      }

      ## Covariance results

      VCVPartitionResults$covariancePartition <-  data.frame("Trait_1" = trait1,
                                                             "Trait_2" = trait2,
                                                             "N" = length(modellingData$dta$animal),

                                                             "Total_coordination" = mean(total_covariance),
                                                             "Total_coordinated_phylogenetic_conservatism" = mean(total_phylogenetic_covariance),
                                                             "Total_coordinated_radiation" = mean(total_non_phylogenetic_covariance),
                                                             "Pvalue_Total_coordination" = total_covariance_pvalue,
                                                             "Pvalue_Total_coordinated_phylogenetic_conservatism" = total_phylogenetic_covariance_pvalue,
                                                             "Pvalue_Total_coordinated_radiation" = total_non_phylogenetic_covariance_pvalue
      )

      VCVPartitionResults$covariancePartitionDistributions <- list("totalCoordination" = total_covariance,
                                                                   "totalCoordinatedPhylogeneticConservatism" = total_phylogenetic_covariance,
                                                                   "Total_coordinated_radiation" = total_non_phylogenetic_covariance)


      if(!is.null(environmental.variable)){
        VCVPartitionResults$covariancePartition <-  cbind(VCVPartitionResults$covariancePartition,
                                                          "Environmental_variables" = paste0(environmental.variable, collapse = ", "),
                                                          "Pure_coordinated_phylogenetic_conservatism" = mean(pure_phylogenetic_covariance),
                                                          "Coordinated_phylogenetic_niche_conservatism" = mean(environmental_phylogenetic_covariance),
                                                          "Pure_environmental_coordination" = mean(pure_environmental_covariance),
                                                          "Residual_coordination" = mean(pure_residual_covariance),

                                                          "Pvalue_Pure_coordinated_phylogenetic_conservatism" = pure_phylogenetic_covariance_pvalue,
                                                          "Pvalue_Coordinated_phylogenetic_niche_conservatism" = environmental_phylogenetic_covariance_pvalue,
                                                          "Pvalue_Pure_environmental_coordination" = pure_environmental_covariance_pvalue,
                                                          "Pvalue_Residual_coordination" = pure_residual_covariance_pvalue
        )

        VCVPartitionResults$covariancePartitionDistributions[["pureCoordinatedPhylogeneticConservatism"]] <- pure_phylogenetic_covariance
        VCVPartitionResults$covariancePartitionDistributions[["coordinatedPhylogeneticNicheConservatism"]] <- environmental_phylogenetic_covariance
        VCVPartitionResults$covariancePartitionDistributions[["pureEnvironmentalCoordination"]] <- pure_environmental_covariance
        VCVPartitionResults$covariancePartitionDistributions[["residualCoordination"]] <- pure_residual_covariance
      }

      VCVPartitionResults$modelPhyloEnv <- mdlPhyloEnv
      VCVPartitionResults$model.diagnostics <- model.diagnostics

      # add to all traits results
      traitsVCVPartitionResults$covarianceResults <- rbind(traitsVCVPartitionResults$covarianceResults,
                                                           VCVPartitionResults$covariancePartition)

      traitsVCVPartitionResults$varianceResults <-  rbind(traitsVCVPartitionResults$varianceResults,
                                                          VCVPartitionResults$variancePartition)

      traitsVCVPartitionResults$models.diagnostics <- rbind(traitsVCVPartitionResults$models.diagnostics,
                                                            VCVPartitionResults$model.diagnostics)
      traitsVCVPartitionResults$individual.models.results[[model]] <- VCVPartitionResults

    } # end evaluation if model already exists in results

  } # end bucle for all traits

  print("Model structure used:")
  print(multi_mdls.str)
  print("Variance results:")
  print(traitsVCVPartitionResults$varianceResults)
  print("Covariance results:")
  print(traitsVCVPartitionResults$covarianceResults)

  # save results

  if(save){
    if(!is.null(environmental.variable)){
      assign(paste0("traitsVCVPartitionResults_", environmental.variable), traitsVCVPartitionResults)
      save(list = paste0("traitsVCVPartitionResults_", environmental.variable), file = results.file)
      print(results.file)
    } else{
      save(list = paste0("traitsVCVPartitionResults"), file = results.file)
      print(results.file)
    }
  }
  return(traitsVCVPartitionResults)

}
