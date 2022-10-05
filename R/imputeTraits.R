#' Impute traits using the evolutionary (or genetic) correlations.
#'
#' @param imputationVariables (character) Names of the variables with NA where imputations will be implemented. If more than one, covariation among vimputation variables is considered.
#' @param dataset (data frame) dataset containing the variable of interest and a column named taxon describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) phylogeny with tip labels contained in dataset$taxon
#' @param correlationsTraitsResults (list) Results from computeCovariancePartition() function or NULL to compute it internally.
#' @param varianceResults (lis)Results from computeVariancePartition() function.
#' @param orderCriterium (character) Name of the correlation to be used as ordering criterium (one of teh columns of correlationsTraitsResults)
#' @param numberOfPhyloCoordinates (number) Number of phylogenetic axis to include
#' @param predictors (character) Names of the variables without NAs used as predictors.
#' @param prodNAs (numeric) Proportion of artificial NAs to be introduced in a given dataset. For evaluation puposes.
#' @param IterationsNumber (integer) Number of iterations of the imputation processs. Results are summarized.
#' @param clustersNumber (integer) Number of clusters to use in parallelization.
#' @param forceRun (logical) If false, previously calculated phylogenetic eigenvectors are used.
#' @param parallelization (logical) If false, no parallelization is conducted.
#' @param numberOfImputationRounds (integer) Number of imputation rounds to implement.
#'
#' @return
#' @export
#'
#' @examples
imputeTraits <- function(imputationVariables, dataset, phylogeny, correlationsTraitsResults = NULL, varianceResults = NULL, orderCriterium = "Total_coordination",
                         numberOfPhyloCoordinates = NULL, predictors = NULL, prodNAs = 0, IterationsNumber = 10, clustersNumber = 2, forceRun = T,
                         parallelization = T, specifications = NULL, numberOfImputationRounds = 3){


  ## only include variables that present a relatively high relationship with the variable to be predicted

  filterVariablesToBeIncludedAsPredictors <- function(imputationVariable, imputedVariables = "", potentialPredictors = NULL, includePhylo = T,
                                                      nPhyloCoords = numberOfPhyloCoordinates){

    # Phylogeny

    if(includePhylo){
      phyloPreds <- character()

      # only for phylogenetically traits
      if(varianceResults[varianceResults$Trait == imputationVariable, "Pvalue_Total_phylogenetic_conservatism"] < 0.05){
        phyloPreds <- c(paste0("Phylo_axis_", 1:nPhyloCoords))
      }
    }

    # Environemntal variables

    if(!is.null(potentialPredictors)){

      correlationsTraitsResults$order <- abs(correlationsTraitsResults[, "Total_coordination"])
      predictors_order <- correlationsTraitsResults %>%
        dplyr::arrange(dplyr::desc(order)) %>%
        dplyr::select(c(Trait_1, Trait_2, Total_coordination, Pvalue_Total_coordination)) %>%
        dplyr::filter(Trait_1 == imputationVariable | Trait_2 == imputationVariable) %>%
        dplyr::filter(Trait_1 %in% potentialPredictors | Trait_2 %in% potentialPredictors) %>%
        dplyr::filter(Pvalue_Total_coordination < 0.05)

      # if(length(predictors_order) < 1){
      # print("not significant effects, selecting those with the highest correlation as predictors (|corr| > Q3)")
      #   threshold <- summary(abs(predictors_order[, "Total_coordination"]))[5]
      #
      #   predictors_order <- predictors_order[predictors_order[, "Total_coordination"] >= threshold | predictors_order[, "Total_coordination"] <= -threshold, ]
      # }

      suitablePredictors <- c(predictors_order$Trait_1, predictors_order$Trait_2)
      suitablePredictors <- unlist(predictors_order[predictors_order %in% potentialPredictors])
      names(suitablePredictors) <- NULL
    } else{
      suitablePredictors <- character()
    }


    # Traits
    suitableTraits_1 <- phyloCov_order[phyloCov_order$Trait_1 == imputationVariable, c("Trait_2", orderCriterium,  paste0("Pvalue_", orderCriterium))] %>%
      dplyr::rename(trait = Trait_2)
    suitableTraits_2 <- phyloCov_order[phyloCov_order$Trait_2 == imputationVariable, c("Trait_1", orderCriterium,  paste0("Pvalue_", orderCriterium))] %>%
      dplyr::rename(trait = Trait_1)
    suitableTraits <- rbind(suitableTraits_1, suitableTraits_2)

    suitableVariables <- character()
    suitableVariables <- suitableTraits[suitableTraits[, paste0("Pvalue_", orderCriterium)] <= 0.5, "trait"]

    # if(length(suitableVariables) < 1){
    #   print("not significant effects, selecting those with the highest correlation as predictors (|corr| > Q3)")
    #   threshold <- summary(abs(suitableTraits[, orderCriterium]))[5]
    #   suitableVariables <- suitableTraits[suitableTraits[, orderCriterium] >= threshold | suitableTraits[, orderCriterium] <= -threshold, "trait"]
    #   suitableVariables <- suitableTraits[suitableTraits[, orderCriterium] >= threshold | suitableTraits[, orderCriterium] <= -threshold, "trait"]
    # }

    imputedSuitableVariables <- suitableVariables[suitableVariables %in% imputedVariables]

    predictiveVariables <- c(imputedSuitableVariables, phyloPreds, suitablePredictors)

    if(length(predictiveVariables) < 1){
      warning("No strong relationships with phylogeny or environment. Using phylogeny to predict.")
      predictiveVariables <- c(paste0("Phylo_axis_", 1:numberOfPhyloCoordinates))
    }

    return(predictiveVariables)
  }

  ### Objects to gather results ####

  predictivePerformance.all <- data.frame()
  ximp.all <- data.frame("taxon" = dataset$taxon)

  predictivePerformance.all2 <- data.frame()
  ximp.all2 <- data.frame("taxon" = dataset$taxon)

  imp.dataset.list <- list()

  imputationResults <- list()

  #### VARIANCE COVARIANCE RESULTS --------------------------------------------- ####

  if(is.null(varianceResults) | is.null(correlationsTraitsResults)){

    if(is.null(specifications)){
      specifications <- defineModelsSpecifications()
    }

    dataset$animal <- dataset$taxon

    vcvResults <- computeVarianceCovariancePartition(traits = c(imputationVariables, predictors), dataset = dataset,
                                                     phylogeny = phylogeny, model.specifications = specifications, force.run = T, save = F,
                                                     showRelativeResults = T, verbose = F)

    if(is.null(varianceResults)){
      varianceResults <- vcvResults$varianceResults
    }

    if(is.null(correlationsTraitsResults)){

      correlationsTraitsResults <- vcvResults$covarianceResults
    }
  }



  #### PHYLOGENETIC PRNCIPAL COMPONENTS ---------------------------------------- ####


  if (!file.exists(paste0(outputs.dir, "phylo_eigenvectors.csv")) | isTRUE(forceRun)) {
    mat <- ape::cophenetic.phylo(phylogeny)
    mat <- as.data.frame(mat)
    pCordA <- ape::pcoa(mat)
    phyloEigenV <- as.data.frame(pCordA$vectors)

    if(is.null(numberOfPhyloCoordinates)){
      numberOfPhyloCoordinates <- max(which(pCordA$values$Relative_eig > 0.01))
    }

    phyloEigenV <- phyloEigenV[dataset$taxon, 1:numberOfPhyloCoordinates] #c(1:50)
    phyloEigenV$taxon <- rownames(phyloEigenV)
    data <- merge(dataset, phyloEigenV, by = "taxon", all.x = T)
    names(data) <- stringr::str_replace_all(names(data), "Axis.", "Phylo_axis_")
    utils::write.csv(data, paste0(outputs.dir, "phylo_eigenvectors.csv"),
                     row.names = F)
    print(paste0(outputs.dir, "phylo_eigenvectors.csv"))
  } else {
    data <- utils::read.csv(paste0(outputs.dir, "phylo_eigenvectors.csv"),
                            header = T)

    if(is.null(numberOfPhyloCoordinates)){
      numberOfPhyloCoordinates <- length(colnames(data)[stringr::str_detect(colnames(data), "Phylo_axis_")])
    }

    print("loading previously calculated phylogenetic eigenvectors")
  }

  if(numberOfPhyloCoordinates > length(colnames(data)[stringr::str_detect(colnames(data), "Phylo_axis_")])){
    stop("numberOfPhyloCoordinates > phylogenetic coordinates explaining > 1% of the variance. Decrease numberOfPhyloCoordinates")
  }

  ### Imputation order according to variance ####

  # From highest phylogenetic variance to lowest

  if(!is.null(varianceResults)){
    imputation.variables_1 <-  varianceResults %>%
      dplyr::arrange(dplyr::desc(abs(Total_phylogenetic_conservatism))) %>%
      dplyr::filter(Trait %in% imputationVariables) %>%
      dplyr::pull(Trait)
  } else {
    imputation.variables_1 <- sample(imputationVariables)
  }

  imputation.variables_1 <- imputation.variables_1[imputation.variables_1 %in% imputationVariables]


  ### Imputation order according to covariance ####

  if(!is.null(orderCriterium)){
    correlationsTraitsResults$order <- abs(correlationsTraitsResults[, orderCriterium])
    phyloCov_order <- correlationsTraitsResults %>%
      dplyr::arrange(dplyr::desc(order)) %>%
      dplyr::select(c(Trait_1, Trait_2, orderCriterium, paste0("Pvalue_", orderCriterium))) %>%
      dplyr::filter(Trait_1 %in% imputationVariables, Trait_2 %in% imputationVariables)

    imputation.variables_2 <- character()
    for (i in 1:length(phyloCov_order[, 1])) {
      impVars <- phyloCov_order[i, ] %>% dplyr::select(Trait_1, Trait_2)
      impVars <- c(impVars$Trait_1, impVars$Trait_2)
      impVars <- impVars[!impVars %in% imputation.variables_2]
      imputation.variables_2 <- c(imputation.variables_2, impVars)
    }

    # Order the first two traits according to their phylogenetic variance
    if(!is.null(varianceResults)){
      firstTwoOrder <-  varianceResults %>%
        dplyr::filter(Trait %in% imputation.variables_2[1:2]) %>%
        dplyr::arrange(dplyr::desc(abs(Total_phylogenetic_conservatism))) %>%
        dplyr::pull(Trait)

      imputation.variables_2[1:2] <- firstTwoOrder
    }
  } else {
    imputation.variables_2 <- sample(imputationVariables)
  }

  imputation.variables_2 <- imputation.variables_2[imputation.variables_2 %in% imputationVariables]

  # Imputation dataset

  imputedVariables <- character() # to store which variables has been imputed already

  ### IMPUTATION 1. Run models to impute variables following the previously determined order using phylogenetic and environmental data ####

  for(imputationVariable in imputation.variables_1){

    if(!is.null(varianceResults)){
      predictiveVariables <- filterVariablesToBeIncludedAsPredictors(imputationVariable = imputationVariable, potentialPredictors = predictors,
                                                                     includePhylo = T, nPhyloCoords = numberOfPhyloCoordinates)
    } else{
      predictiveVariables <- c(predictors, c(paste0("Phylo_axis_", 1:numberOfPhyloCoordinates)))
    }

    modelName <- paste("(", paste(imputationVariable, collapse = " <-> "), ") <- ", paste(predictiveVariables, collapse = ", "))

    imp.dataset <- data %>% dplyr::select(c(imputationVariable, predictiveVariables))
    xTrue <- imp.dataset

    for (n in 1:IterationsNumber) {

      imp.dataset[, imputationVariable] <- missForest::prodNA(as.data.frame(xTrue[, imputationVariable]), prodNAs)

      imp.dataset.list[[imputationVariable]][[n]] <- imp.dataset

      cl <- parallel::makeCluster(clustersNumber)
      doParallel::registerDoParallel(cl)
      if (prodNAs != 0) {
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
                                            NRMSE = rfImp.res$error[1:length(imputationVariable)],
                                            R2 = r2.var,
                                            Model = modelName)
      } else {

        if(all(!is.na(imp.dataset))){
          stop("No missing values in dataset")
        }

        rfImp.res <- randomForestImpute(xmis = as.matrix(imp.dataset),
                                        maxiter = 50, ntree = 1000, parallelize = parallelization)

        predictivePerformance <- data.frame(Variable = imputationVariable,
                                            N = length(imp.dataset[, 1]),
                                            N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, imputationVariable])), 1]),
                                            N_NA = length(imp.dataset[which(is.na(imp.dataset[, imputationVariable])), 1]),
                                            NRMSE = rfImp.res$OOBerror[1:length(imputationVariable)],
                                            Model = modelName)
      }
      parallel::stopCluster(cl)

      row.names(predictivePerformance) <- NULL
      print(predictivePerformance)

      ## Results of the first round of imputation (gap filling) per iteration

      ximp <- as.data.frame(rfImp.res$ximp)
      ximp$taxon <- dataset$taxon
      predictivePerformance.all <- rbind(predictivePerformance.all, predictivePerformance)
      ximp.all[, imputationVariable] <- ximp[, imputationVariable]
    }

  }

  ### Results aggregation (for all iterations)
  imputationResults[["round1"]]$modelFormula <- modelName
  imputationResults[["round1"]]$ximp <- stats::aggregate(ximp.all[, -which(names(ximp.all) == "taxon")], by = list(ximp.all$taxon), FUN = mean) %>% dplyr::rename(taxon = Group.1)
  imputationResults[["round1"]]$predictivePerformance <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"))],
                                                                          by = list(predictivePerformance.all$Variable), FUN = meanOrMode) %>%
    dplyr::rename(Variable = Group.1)

  if (IterationsNumber > 1) {
    imputationResults[["round1"]]$ximp_sd <- stats::aggregate(ximp.all[, -which(names(ximp.all) == "taxon")], by = list(ximp.all$taxon), FUN = stats::sd) %>%
      dplyr::rename(taxon = Group.1)

    imputationResults[["round1"]]$predictivePerformance_sd <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"), which(names(predictivePerformance.all) == "Model"))],
                                                                               by = list(predictivePerformance.all$Variable), FUN = stats::sd) %>% dplyr::rename(Variable = Group.1)
    imputationResults[["round1"]]$ximp_all_iterations <- ximp.all
    imputationResults[["round1"]]$predictivePerformance_all_iterations <- predictivePerformance.all
  }


  ### IMPUTATION ROUNDS 2:N. Run models to impute variables following the previously determined order also including imputed traits  ####

  for(nround in 2:numberOfImputationRounds){
    ximp <- as.data.frame(imputationResults[[paste0("round", nround - 1)]]$ximp)
    ximp <- merge(ximp, data[, -which(names(data) %in% imputationVariables)], by = "taxon")

    for(imputationVariable in imputation.variables_2){

      if(!is.null(correlationsTraitsResults)){
        predictiveVariables <- filterVariablesToBeIncludedAsPredictors(imputationVariable = imputationVariable, imputedVariables = imputation.variables_2,
                                                                       potentialPredictors = predictors, includePhylo = T, nPhyloCoords = numberOfPhyloCoordinates)

      } else{
        predictiveVariables <- c(sample(imputation.variables_2), predictors, c(paste0("Phylo_axis_", 1:numberOfPhyloCoordinates)))
      }

      modelName <- paste("(", paste(imputationVariable, collapse = " <-> "), ") <- ", paste(predictiveVariables, collapse = ", "))

      imp.dataset <- ximp %>% dplyr::select(c(imputationVariable, predictiveVariables))

      for(n in 1:IterationsNumber) {

        imp.dataset[, imputationVariable] <- imp.dataset.list[[imputationVariable]][[n]][, imputationVariable] # it already have NAs generated in the previous step (if needed)

        xTrue <- data %>% dplyr::select(c(imputationVariable, predictiveVariables))

        cl <- parallel::makeCluster(clustersNumber)
        doParallel::registerDoParallel(cl)
        if (prodNAs != 0) {
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

          predictivePerformance <- data.frame(Variable = imputationVariable,
                                              N = length(ximp[, 1]),
                                              N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, imputationVariable])), 1]),
                                              N_NA = length(imp.dataset[which(is.na(imp.dataset[, imputationVariable])), 1]),
                                              NRMSE = rfImp.res$error[1:length(imputationVariable)],
                                              R2 = r2.var,
                                              Model = modelName)
        } else {
          rfImp.res <- randomForestImpute(xmis = as.matrix(imp.dataset),
                                          maxiter = 50, ntree = 1000, parallelize = parallelization)

          predictivePerformance <- data.frame(Variable = imputationVariable,
                                              N = length(imp.dataset[, 1]),
                                              N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, imputationVariable])), 1]),
                                              N_NA = length(imp.dataset[which(is.na(imp.dataset[, imputationVariable])), 1]),
                                              NRMSE = rfImp.res$OOBerror[1:length(imputationVariable)],
                                              Model = modelName)
        }
        parallel::stopCluster(cl)

        ## Results of the second round of imputation (gap filling) per iteration

        row.names(predictivePerformance) <- NULL
        print(predictivePerformance)
        imp.dataset[, imputationVariable] <- as.data.frame(rfImp.res$ximp[, imputationVariable])

        predictivePerformance.all <- rbind(predictivePerformance.all, predictivePerformance)
        ximp.all[, imputationVariable] <- imp.dataset[, imputationVariable]
      }
    }

    ### Results aggregation (for all iterations)
    imputationResults[[paste0("round", nround)]]$modelFormula <- modelName

    imputationResults[[paste0("round", nround)]]$ximp <- stats::aggregate(ximp.all[, -which(names(ximp.all) == "taxon")], by = list(ximp.all$taxon), FUN = mean) %>%
      dplyr::rename(taxon = Group.1)
    imputationResults[[paste0("round", nround)]]$predictivePerformance <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"))],
                                                                                           by = list(predictivePerformance.all$Variable), FUN = meanOrMode) %>%
      dplyr::rename(Variable = Group.1)

    if (IterationsNumber > 1) {
      imputationResults[[paste0("round", nround)]]$ximp_sd <- stats::aggregate(ximp.all[, -which(names(ximp.all) == "taxon")], by = list(ximp.all$taxon), FUN = stats::sd) %>%
        dplyr::rename(taxon = Group.1)

      imputationResults[[paste0("round", nround)]]$predictivePerformance_sd <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"), which(names(predictivePerformance.all) == "Model"))],
                                                                                                by = list(predictivePerformance.all$Variable), FUN = stats::sd) %>% dplyr::rename(Variable = Group.1) %>%
        dplyr::mutate(Model = modelName)
      imputationResults[[paste0("round", nround)]]$ximp_all_iterations <- ximp.all
      imputationResults[[paste0("round", nround)]]$predictivePerformance_all_iterations <- predictivePerformance.all
    }
  }

  return(imputationResults)
}
