#' Impute traits using phylogenetic, environmental and traits relationships
#'
#' Predict missing values using phylogeny, environmental and trait data as predictors in random forests when it shows a relationship with the variable
#' to be imputed. The methodology implements 3 rounds in which the actualized predicted values are used to predict the rest of the variables to impute.
#'
#' @param variables_to_impute (*character*). Names of the variables with NAs where imputations will be implemented.
#' If more than one, covariation among imputation variables is considered to decide whether to use them as predictors or not.
#' @param dataset (*data frame*). Data frame containing the variable of interest with missing values and a column describing terminal taxa of phylogeny (e.g., species).
#' @param terminal_taxon (*character*). Terminal taxon as named in the dataset (e.g., species).
#' @param phylogeny (*phylo*). Phylogeny with tip labels contained in dataset.
#' @param correlation_results (*list*). Correlation results from computeVarianceCovariance() function or NULL to compute it internally for the variables to be imputed using the data and phylogeny provided (default).
#' @param variance_results (*list*). Variance results from computeVarianceCovariance() function.
#' @param number_of_phylo_axis (*integer*). Number of phylogenetic axis to include as predictors.
#' @param predictors (*character*). Names of the variables without NAs that will be considered as potential. These varaibles need to be included in the dataset and need to be complete (no NAs).
#' @param proportion_NAs (*numeric*). Between 0 and 1. Proportion of artificial NAs to be introduced in a given dataset. For evaluation puposes. Default set to 0, so no extra NAs are produced.
#' @param number_iterations (*integer*). Number of iterations of the imputation process. Results reported are summarized as mean and standard deviation for each variable to be imputed.
#' @param number_clusters (*integer*). Number of clusters to use in parallelization. If set to one, no paralelization is performed (not recommended). Default is set to 2.
#' @param model_specifications (*list*). Mcmcglmm models specifications as specified by the defineModelsSpecification of this package. If not defined,
#' the  default of the defineModelsSpecification function is used.
#' @param save (*logical*) If false, results are not saved in the outputs folder.
#' @param force_run (*logical*) If false, previously calculated phylogenetic axis (saved internally) are used. This can be useful when using big phylogenies, as the calculation
#' of the phylogenetic axis can take some time.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate example data
#' simulated_traits.data <- simulateDataSet()
#'
#' # Impute missing values (created within the function)
#' imputed.data <- imputetraits(
#' variables_to_impute = c("phylo_G1_trait1", "phylo_G1_trait2"),
#' dataset = simulated_traits.data$data,
#' phylogeny = simulated_traits.data$phylogeny,
#' proportion_NAs = 0.2
#' )
#' }
imputeTraits <- function(variables_to_impute = NULL,
                         dataset = NULL,
                         terminal_taxon = NULL,
                         phylogeny = NULL,
                         correlation_results = NULL,
                         variance_results = NULL,
                         number_of_phylo_axis = NULL,
                         predictors = NULL,
                         proportion_NAs = 0,
                         number_iterations = 10,
                         number_clusters = 2,
                         model_specifications = NULL,
                         save = F,
                         force_run = T
                         ){

  # Arguments
  if(is.null(variables_to_impute)){
    stop("Specify variables_to_impute argument")
  }

  if(is.null(dataset)){
    stop("Specify dataset argument")
  }

  if(is.null(phylogeny)){
    stop("Specify phylogeny argument")
  }

  if(is.null(terminal_taxon)){
    stop("Specify terminal_taxon argument")
  }

  # Check whether initializeTrevol has been run if user wants to automatically save results
  if(isTRUE(save)){
    if(isFALSE(exists("outputs.dir"))){
      stop("If you set save = T you first need to run initializeTrEvol to create the folder and subfolder structure where results are saved.")
    }
  }


  # Create taxon column
  dataset$taxon <- dataset[, terminal_taxon]

  # Set correlation type to total correlation
  correlation_type <- "total_correlation"

  # Set number of imputation rounds to 3
  numberOfImputationRounds <- 3

  # Set parallelization argument internally
  if(number_clusters > 1){
    parallelization <- T
  } else{
    parallelization <- F
  }

  ## only include variables that present a relatively high relationship with the variable to be predicted

  filterVariablesToBeIncludedAsPredictors <- function(imputationVariable,
                                                      imputedVariables = "",
                                                      potentialPredictors = NULL,
                                                      includePhylo = T,
                                                      nPhyloCoords = number_of_phylo_axis){

    # Phylogeny

    if(includePhylo){
      phyloPreds <- character()

      # only for phylogenetically traits
      if(variance_results[variance_results$trait == imputationVariable, "CI_significance_phylogenetic_variance"] == "yes"){
        phyloPreds <- c(paste0("Phylo_axis_", 1:nPhyloCoords))
      }
    }

    # Environemntal variables

    if(!is.null(potentialPredictors)){

      correlation_results$order <- abs(correlation_results[, "total_correlation"])
      predictors_order <- correlation_results %>%
        dplyr::arrange(dplyr::desc(order)) %>%
        dplyr::select(c(trait_1, trait_2, total_correlation, CI_significance_total_correlation)) %>%
        dplyr::filter(trait_1 == imputationVariable | trait_2 == imputationVariable) %>%
        dplyr::filter(trait_1 %in% potentialPredictors | trait_2 %in% potentialPredictors) %>%
        dplyr::filter(CI_significance_total_correlation == "yes")

      # if(length(predictors_order) < 1){
      # print("not significant effects, selecting those with the highest correlation as predictors (|corr| > Q3)")
      #   threshold <- summary(abs(predictors_order[, "total_correlation"]))[5]
      #
      #   predictors_order <- predictors_order[predictors_order[, "total_correlation"] >= threshold | predictors_order[, "total_correlation"] <= -threshold, ]
      # }

      suitablePredictors <- c(predictors_order$trait_1, predictors_order$trait_2)
      suitablePredictors <- unlist(predictors_order[predictors_order %in% potentialPredictors])
      names(suitablePredictors) <- NULL
    } else{
      suitablePredictors <- character()
    }


    # traits
    suitabletraits_1 <- phyloCov_order[phyloCov_order$trait_1 == imputationVariable, c("trait_2", correlation_type,  paste0("CI_significance_", correlation_type))] %>%
      dplyr::rename(trait = trait_2)
    suitabletraits_2 <- phyloCov_order[phyloCov_order$trait_2 == imputationVariable, c("trait_1", correlation_type,  paste0("CI_significance_", correlation_type))] %>%
      dplyr::rename(trait = trait_1)
    suitabletraits <- rbind(suitabletraits_1, suitabletraits_2)

    suitableVariables <- character()
    suitableVariables <- suitabletraits[suitabletraits[, paste0("CI_significance_", correlation_type)] == "yes", "trait"]

    # if(length(suitableVariables) < 1){
    #   print("not significant effects, selecting those with the highest correlation as predictors (|corr| > Q3)")
    #   threshold <- summary(abs(suitabletraits[, correlation_type]))[5]
    #   suitableVariables <- suitabletraits[suitabletraits[, correlation_type] >= threshold | suitabletraits[, correlation_type] <= -threshold, "trait"]
    #   suitableVariables <- suitabletraits[suitabletraits[, correlation_type] >= threshold | suitabletraits[, correlation_type] <= -threshold, "trait"]
    # }

    imputedSuitableVariables <- suitableVariables[suitableVariables %in% imputedVariables]

    predictiveVariables <- c(imputedSuitableVariables, phyloPreds, suitablePredictors)

    if(length(predictiveVariables) < 1){
      warning("No strong relationships with phylogeny or environment. Using phylogeny to predict.")
      predictiveVariables <- c(paste0("Phylo_axis_", 1:number_of_phylo_axis))
    }

    return(predictiveVariables)
  }

  ### Objects to gather results ####

  predictivePerformance.all <- data.frame()

  ## Create an ID column that will be used to aggregate results by row later

  ximp.all <- data.frame("taxon" = rep(dataset$taxon, number_iterations),
                         "ID" = rep(1:dim(dataset)[1], number_iterations))

  predictivePerformance.all_copy <- predictivePerformance.all
  ximp.all_copy <- ximp.all

  imp.dataset.list <- list()
  imputationResults <- list()


  #### VARIANCE COVARIANCE RESULTS --------------------------------------------- ####

  if(is.null(variance_results) | is.null(correlation_results)){

    if(is.null(model_specifications)){
      model_specifications <- defineModelsSpecifications()
    }

    dataset$animal <- dataset$taxon

    vcvResults <- computeVarianceCovariancePartition(traits = c(variables_to_impute, predictors),
                                                     dataset = dataset,
                                                     phylogeny = phylogeny,
                                                     model_specifications = model_specifications,
                                                     terminal_taxa = terminal_taxon)

    if(is.null(variance_results)){
      variance_results <- vcvResults$varianceResults
    }

    if(is.null(correlation_results)){

      correlation_results <- vcvResults$covarianceResults
    }
  }


  #### PHYLOGENETIC PRNCIPAL COMPONENTS ---------------------------------------- ####

  if(isTRUE(force_run)) {
    mat <- ape::cophenetic.phylo(phylogeny)
    mat <- as.data.frame(mat)
    pCordA <- ape::pcoa(mat)
    phyloEigenV <- as.data.frame(pCordA$vectors)

    if(is.null(number_of_phylo_axis)){
      number_of_phylo_axis <- max(which(pCordA$values$Relative_eig > 0.01))
    }

    phyloEigenV <- phyloEigenV[dataset$taxon, 1:number_of_phylo_axis] #c(1:50)
    phyloEigenV$taxon <- rownames(phyloEigenV)
    data <- merge(dataset, phyloEigenV, by = "taxon", all.x = T)
    names(data) <- stringr::str_replace_all(names(data), "Axis.", "Phylo_axis_")

    if(save){
      utils::write.csv(data, paste0(outputs.dir, "phylo_eigenvectors.csv"),
                       row.names = F)
      message(paste0(outputs.dir, "phylo_eigenvectors.csv"))
    }
  } else {
    data <- utils::read.csv(paste0(outputs.dir, "phylo_eigenvectors.csv"),
                            header = T)

    if(is.null(number_of_phylo_axis)){
      number_of_phylo_axis <- length(colnames(data)[stringr::str_detect(colnames(data), "Phylo_axis_")])
    }

    message("loading previously calculated phylogenetic eigenvectors")
  }

  if(number_of_phylo_axis > length(colnames(data)[stringr::str_detect(colnames(data), "Phylo_axis_")])){
    stop("number_of_phylo_axis > phylogenetic coordinates explaining > 1% of the variance. Decrease number_of_phylo_axis")
  }

  ### Imputation order according to variance ####

  # From highest phylogenetic variance to lowest

  if(!is.null(variance_results)){
    imputation.variables_1 <-  variance_results %>%
      dplyr::arrange(dplyr::desc(abs(phylogenetic_variance))) %>%
      dplyr::filter(trait %in% variables_to_impute) %>%
      dplyr::pull(trait)
  } else {
    imputation.variables_1 <- sample(variables_to_impute)
  }

  imputation.variables_1 <- imputation.variables_1[imputation.variables_1 %in% variables_to_impute]


  ### Imputation order according to covariance ####

  if(!is.null(correlation_type)){
    correlation_results$order <- abs(correlation_results[, correlation_type])

    if(length(correlation_results[, correlation_type]) > 1){
      phyloCov_order <- correlation_results %>%
        dplyr::arrange(dplyr::desc(order)) %>%
        dplyr::select(c(trait_1, trait_2, all_of(correlation_type), paste0("CI_significance_", correlation_type))) %>%
        dplyr::filter(trait_1 %in% variables_to_impute, trait_2 %in% variables_to_impute)

      imputation.variables_2 <- character()
      for (i in 1:length(phyloCov_order[, 1])) {
        impVars <- phyloCov_order[i, ] %>% dplyr::select(trait_1, trait_2)
        impVars <- c(impVars$trait_1, impVars$trait_2)
        impVars <- impVars[!impVars %in% imputation.variables_2]
        imputation.variables_2 <- c(imputation.variables_2, impVars)
      }
    } else{
      phyloCov_order <- correlation_results

      imputation.variables_2 <- c(correlation_results$trait_1, correlation_results$trait_2)
    }



    # Order the first two traits according to their phylogenetic variance
    if(!is.null(variance_results)){
      firstTwoOrder <-  variance_results %>%
        dplyr::filter(trait %in% imputation.variables_2[1:2]) %>%
        dplyr::arrange(dplyr::desc(abs(phylogenetic_variance))) %>%
        dplyr::pull(trait)

      imputation.variables_2[1:2] <- firstTwoOrder
    }
  } else {
    imputation.variables_2 <- sample(variables_to_impute)
  }

  imputation.variables_2 <- imputation.variables_2[imputation.variables_2 %in% variables_to_impute]

  # Imputation dataset

  imputedVariables <- character() # to store which variables has been imputed already


  ### IMPUTATION 1. Run models to impute variables following the previously determined order using phylogenetic and environmental data ####

  for(imputationVariable in imputation.variables_1){

    if(!is.null(variance_results)){
      predictiveVariables <- filterVariablesToBeIncludedAsPredictors(imputationVariable = imputationVariable, potentialPredictors = predictors,
                                                                     includePhylo = T, nPhyloCoords = number_of_phylo_axis)
    } else{
      predictiveVariables <- c(predictors, c(paste0("Phylo_axis_", 1:number_of_phylo_axis)))
    }

    modelName <- paste("(", paste(imputationVariable, collapse = " <-> "), ") <- ", paste(predictiveVariables, collapse = ", "))

    imp.dataset <- data %>% dplyr::select(c(all_of(imputationVariable), all_of(predictiveVariables)))
    xTrue <- imp.dataset

    # Data frames to get imputation variable results
    ximp <- data.frame()
    predictivePerformance <- data.frame()

    for (n in 1:number_iterations) {

      # Produce NAs if proportion_NAs is not zero
      imp.dataset[, imputationVariable] <- missForest::prodNA(as.data.frame(xTrue[, imputationVariable]), proportion_NAs)

      # Save variable to be imputed with NAs generated to calculate performance in next iterations round
      imp.dataset.list[[imputationVariable]][[n]] <- imp.dataset

      cl <- parallel::makeCluster(number_clusters)
      doParallel::registerDoParallel(cl)
      if (proportion_NAs != 0) {
        rfImp.res <- randomForestImpute(xmis = as.matrix(imp.dataset),
                                        maxiter = 50, ntree = 100, parallelize = parallelization,
                                        xtrue = as.matrix(xTrue))

        # R2 calculation
        r2.var <- numeric()
        matObs <- as.matrix(xTrue[, imputationVariable])
        matNA <- as.matrix(imp.dataset[, imputationVariable])

        observed <- matObs[which(!is.na(matObs) & is.na(matNA))]

        predicted <- as.matrix(rfImp.res$ximp[, imputationVariable])
        predicted <- predicted[which(!is.na(matObs) & is.na(matNA))]

        r2.var[imputationVariable] <- cor(observed, predicted)^2


        predictivePerformance <- data.frame(Variable = imputationVariable,
                                            N = length(imp.dataset[, 1]),
                                            N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, imputationVariable])), 1]),
                                            N_NA = length(imp.dataset[which(is.na(imp.dataset[, imputationVariable])), 1]),
                                            NRMSE = rfImp.res[[2]][1:length(imputationVariable)],
                                            R2 = r2.var,
                                            Model = modelName)
      } else {

        if(all(!is.na(imp.dataset))){
          stop("No missing values in dataset")
        }

        rfImp.res <- randomForestImpute(xmis = as.matrix(imp.dataset),
                                        maxiter = 50, ntree = 1000, parallelize = parallelization)

        predictivePerformance_n <- data.frame(Variable = imputationVariable,
                                            N = length(imp.dataset[, 1]),
                                            N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, imputationVariable])), 1]),
                                            N_NA = length(imp.dataset[which(is.na(imp.dataset[, imputationVariable])), 1]),
                                            NRMSE = rfImp.res[[2]][1:length(imputationVariable)],
                                            Model = modelName)
      }
      parallel::stopCluster(cl)

      row.names(predictivePerformance_n) <- NULL
      # print(predictivePerformance_n)

      ## Results of the first round of imputation (gap filling) per iteration

      ximp <- rbind(ximp, as.data.frame(rfImp.res$ximp))
      predictivePerformance <- rbind(predictivePerformance, predictivePerformance_n)
    }

    # Add imputed variable results to the general output
    ximp.all[, imputationVariable] <- ximp[, imputationVariable]
    predictivePerformance.all <- rbind(predictivePerformance.all, predictivePerformance)
  }

  ### Results aggregation (for all iterations)
  imputationResults[["round1"]]$modelFormula <- modelName

  ximp_1 <- stats::aggregate(ximp.all[, variables_to_impute, drop = F], by = list(ximp.all$ID), FUN = mean) %>% dplyr::rename(ID = Group.1)

  imputationResults[["round1"]]$ximp <- merge(ximp.all[1:dim(data)[1], c("taxon", "ID")], ximp_1, by = "ID", all.x = T) %>%
    dplyr::select(taxon, all_of(variables_to_impute)) %>%
    dplyr::arrange(taxon)

  imputationResults[["round1"]]$predictivePerformance <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"))],
                                                                          by = list(predictivePerformance.all$Variable), FUN = meanOrMode) %>%
    dplyr::rename(Variable = Group.1) %>%
    dplyr::mutate(Model = modelName)

  if (number_iterations > 1) {
    ximp_1_sd <- stats::aggregate(ximp.all[, variables_to_impute, drop = F], by = list(ximp.all$ID), FUN = sd) %>% dplyr::rename(ID = Group.1)

    imputationResults[["round1"]]$ximp_sd <- merge(ximp.all[1:dim(data)[1], c("taxon", "ID")], ximp_1_sd, by = "ID", all.x = T) %>%
      dplyr::select(taxon, all_of(variables_to_impute)) %>% dplyr::arrange(taxon)

    imputationResults[["round1"]]$predictivePerformance_sd <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"), which(names(predictivePerformance.all) == "Model"))],
                                                                               by = list(predictivePerformance.all$Variable), FUN = stats::sd) %>%
      dplyr::rename(Variable = Group.1) %>%
      dplyr::mutate(Model = modelName)
    imputationResults[["round1"]]$ximp_all_iterations <- ximp.all
    imputationResults[["round1"]]$predictivePerformance_all_iterations <- predictivePerformance.all
  }


  ### IMPUTATION ROUNDS 2:N. Run models to impute variables following the previously determined order also including imputed traits  ####

  for(nround in 2:numberOfImputationRounds){

    ximp.all2 <- ximp.all_copy
    predictivePerformance.all2 <- predictivePerformance.all_copy

    for(imputationVariable in imputation.variables_2){

      if(!is.null(correlation_results)){
        predictiveVariables <- filterVariablesToBeIncludedAsPredictors(imputationVariable = imputationVariable, imputedVariables = imputation.variables_2,
                                                                       potentialPredictors = predictors, includePhylo = T, nPhyloCoords = number_of_phylo_axis)

      } else{
        predictiveVariables <- c(sample(imputation.variables_2), predictors, c(paste0("Phylo_axis_", 1:number_of_phylo_axis)))
      }

      ximp <- as.data.frame(imputationResults[[paste0("round", nround - 1)]]$ximp)
      ximp <- cbind(ximp, data[, -which(names(data) %in% variables_to_impute)])

      modelName <- paste("(", paste(imputationVariable, collapse = " <-> "), ") <- ", paste(predictiveVariables, collapse = ", "))

      imp.dataset <- ximp %>% dplyr::select(c(all_of(imputationVariable), all_of(predictiveVariables)))


      # Data frames to get imputation variable results
      ximp <- data.frame()
      predictivePerformance <- data.frame()

      for(n in 1:number_iterations) {

        # select the variable to be imputed from the original dataset
        imp.dataset[, imputationVariable] <- imp.dataset.list[[imputationVariable]][[n]][, imputationVariable] # it already have NAs generated in the previous step (if needed)

        xTrue <- data %>% dplyr::select(c(all_of(imputationVariable), all_of(predictiveVariables)))

        cl <- parallel::makeCluster(number_clusters)
        doParallel::registerDoParallel(cl)
        if (proportion_NAs != 0) {
          rfImp.res <- randomForestImpute(xmis = as.matrix(imp.dataset),
                                          maxiter = 50, ntree = 1000, parallelize = parallelization,
                                          xtrue = as.matrix(xTrue))


          # R2 calculation
          matObs <- as.matrix(xTrue[, imputationVariable])
          matNA <- as.matrix(imp.dataset[, imputationVariable])

          observed <- matObs[which(!is.na(matObs) & is.na(matNA))]
          predicted <- as.matrix(rfImp.res$ximp[, imputationVariable])
          predicted <- predicted[which(!is.na(matObs) & is.na(matNA))]

          r2.var <- cor(observed, predicted)^2

          predictivePerformance_n <- data.frame(Variable = imputationVariable,
                                              N = length(ximp[, 1]),
                                              N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, imputationVariable])), 1]),
                                              N_NA = length(imp.dataset[which(is.na(imp.dataset[, imputationVariable])), 1]),
                                              NRMSE = rfImp.res[[2]][1:length(imputationVariable)],
                                              R2 = r2.var,
                                              Model = modelName)
        } else {
          rfImp.res <- randomForestImpute(xmis = as.matrix(imp.dataset),
                                          maxiter = 50, ntree = 1000, parallelize = parallelization)

          predictivePerformance_n <- data.frame(Variable = imputationVariable,
                                              N = length(imp.dataset[, 1]),
                                              N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, imputationVariable])), 1]),
                                              N_NA = length(imp.dataset[which(is.na(imp.dataset[, imputationVariable])), 1]),
                                              NRMSE = rfImp.res[[2]][1:length(imputationVariable)],
                                              Model = modelName)
        }
        parallel::stopCluster(cl)

        ## Results of the second round of imputation (gap filling) per iteration

        row.names(predictivePerformance_n) <- NULL
        # print(predictivePerformance_n)

        # add imputed variable to the imputatin dataset (may be considered as predictor)

        imp.dataset[, imputationVariable] <- as.data.frame(rfImp.res$ximp[, imputationVariable])

        ximp <- rbind(ximp, as.data.frame(rfImp.res$ximp))
        predictivePerformance <- rbind(predictivePerformance, predictivePerformance_n)

      }

      # Add imputed variable results to the general output
      ximp.all2[, imputationVariable] <- ximp[, imputationVariable]
      predictivePerformance.all2 <- rbind(predictivePerformance.all2, predictivePerformance)

    }

    ### Results aggregation (for all iterations)
    imputationResults[[paste0("round", nround)]]$modelFormula <- modelName


    ximp_1 <- stats::aggregate(ximp.all2[, variables_to_impute, drop = F], by = list(ximp.all2$ID), FUN = mean) %>% dplyr::rename(ID = Group.1)

    imputationResults[[paste0("round", nround)]]$ximp <- merge(ximp.all2[1:dim(data)[1], c("taxon", "ID")], ximp_1, by = "ID", all.x = T) %>%
      dplyr::select(taxon, all_of(variables_to_impute)) %>%
      dplyr::arrange(taxon)

    imputationResults[[paste0("round", nround)]]$predictivePerformance <- stats::aggregate(predictivePerformance.all2[, -c(which(names(predictivePerformance.all2) == "Variable"))],
                                                                                           by = list(predictivePerformance.all2$Variable), FUN = meanOrMode) %>%
      dplyr::rename(Variable = Group.1)

    if (number_iterations > 1) {
      ximp_1_sd <- stats::aggregate(ximp.all2[, variables_to_impute, drop = F], by = list(ximp.all2$ID), FUN = sd) %>% dplyr::rename(ID = Group.1)

      imputationResults[[paste0("round", nround)]]$ximp_sd <- merge(ximp.all2[1:dim(data)[1], c("taxon", "ID")], ximp_1_sd, by = "ID", all.x = T) %>%
        dplyr::select(taxon, all_of(variables_to_impute)) %>% dplyr::arrange(taxon)

      imputationResults[[paste0("round", nround)]]$predictivePerformance_sd <- stats::aggregate(predictivePerformance.all2[, -c(which(names(predictivePerformance.all2) == "Variable"), which(names(predictivePerformance.all2) == "Model"))],
                                                                                                by = list(predictivePerformance.all2$Variable), FUN = stats::sd) %>% dplyr::rename(Variable = Group.1) %>%
        dplyr::mutate(Model = modelName)
      imputationResults[[paste0("round", nround)]]$ximp_all_iterations <- ximp.all2
      imputationResults[[paste0("round", nround)]]$predictivePerformance_all_iterations <- predictivePerformance.all2
    }
  }

  # Rename terminal taxon

  names(imputationResults[["round1"]]$ximp)[1] <- terminal_taxon
  names(imputationResults[["round2"]]$ximp)[1] <- terminal_taxon
  names(imputationResults[["round3"]]$ximp)[1] <- terminal_taxon

  return(imputationResults)
}
