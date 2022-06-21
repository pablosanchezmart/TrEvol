#' Impute traits using the evolutionary (or genetic) correlations.
#'
#' @param imputationVariables (character) Names of the variables with NA where imputations will be implemented. If more than one, covariation among vimputation variables is considered.
#' @param predictors (character) Names of the variables without NAs used as predictors.
#' @param dataset (data frame) dataset containing the variable of interest and a column named taxon describing terminal taxa of phylogeny.
#' @param phylogeny (phylo) phylogeny with tip labels contained in dataset$taxon
#' @param correlationsTraitsResults (list) Results from correlationsTraits() function.
#' @param prodNAs (numeric) Proportion of artificial NAs to be introduced in a given dataset. For evaluation puposes.
#' @param IterationsNumber (integer) Number of iterations of the imputation processs. Results are summarized.
#' @param clustersNumber (integer) Number of clusters to use in parallelization.
#' @param forceRun (logical) If false, previously calculated phylogenetic eigenvectors are used.
#'
#' @return
#' @export
#'
#' @examples
imputeTraits <- function(dataset, phylogeny, correlationsTraitsResults, varianceResults = NULL, imputationVariables, orderCriterium = NULL, numberOfPhyloCoordinates = 5,
                             predictors = NULL, prodNAs = 0, IterationsNumber = 10, clustersNumber = 2, forceRun = T){

  ### Objects to gather results ####

  predictivePerformance.all <- data.frame()
  ximp.all <- data.frame()

  predictivePerformance.all2 <- data.frame()
  ximp.all2 <- data.frame()

  imputationResults <- list()

  modelName <- paste("(", paste(imputationVariables, collapse = " <-> "), ") <- ", paste(predictors, collapse = ", "), paste("numberPhyloCoord_", numberOfPhyloCoordinates))

  #### PHYLOGENETIC PRNCIPAL COMPONENTS ------------------------------------------ ####

  ### Phylogenetic principal coordinates ####

  if (!file.exists(paste0(outputs.dir, "phylo_eigenvectors.csv")) | isTRUE(forceRun)) {
    mat <- ape::cophenetic.phylo(phylogeny)
    mat <- as.data.frame(mat)
    pCordA <- ape::pcoa(mat)
    phyloEigenV <- as.data.frame(pCordA$vectors)
    phyloEigenV <- phyloEigenV[dataset$taxon, c(1:50)]
    phyloEigenV$taxon <- rownames(phyloEigenV)
    imp.dataset <- merge(dataset, phyloEigenV, by = "taxon", all.x = T)
    names(imp.dataset) <- stringr::str_replace_all(names(imp.dataset), "Axis.", "Phylo_axis_")
    utils::write.csv(imp.dataset, paste0(outputs.dir, "phylo_eigenvectors.csv"),
                     row.names = F)
    print(paste0(outputs.dir, "phylo_eigenvectors.csv"))
  } else {
    imp.dataset <- utils::read.csv(paste0(outputs.dir, "phylo_eigenvectors.csv"),
                                   header = T)
    print("loading previously calculated phylogenetic eigenvectors")
  }

  ### Imputation order ####

  # Order traits according to the covariance structure

  if(!is.null(orderCriterium)){
  correlationsTraitsResults$order <- abs(correlationsTraitsResults[, orderCriterium])
  phylo_order <- correlationsTraitsResults %>%
    dplyr::arrange(dplyr::desc(order)) %>%
    dplyr::select(c(Trait_1, Trait_2, orderCriterium))

  imputation.variables <- character()
  for (i in 1:length(phylo_order[, 1])) {
    impVars <- phylo_order[i, ] %>% dplyr::select(Trait_1, Trait_2)
    impVars <- c(impVars$Trait_1, impVars$Trait_2)
    impVars <- impVars[!impVars %in% imputation.variables]
    imputation.variables <- c(imputation.variables, impVars)
  }
  } else {
    imputation.variables <- imputationVariables
  }

  # Order the first two traits according to their phylogenetic variance

  firstTwoOrder <-  varianceResults %>%
    dplyr::filter(Trait %in% imputation.variables[1:2]) %>%
    dplyr::arrange(dplyr::desc(abs(Total_phylogenetic_conservatism))) %>%
    dplyr::pull(Trait)

  imputation.variables[1:2] <- firstTwoOrder

  imputationVariables <- imputation.variables[imputation.variables %in% imputationVariables]

  # Imputation dataset

  predictors <- c(predictors, c(paste0("Phylo_axis_", 1:numberOfPhyloCoordinates)))

  imp.dataset <- imp.dataset %>% dplyr::select(c(imputationVariables, predictors))
  xTrue <- imp.dataset

  ### IMPUTATION 1. Run models to impute variables following the previously determined order ####

  for (n in 1:IterationsNumber) {
    # if order criterium is null, randomize order in each iteration
    if(is.null(orderCriterium)){
      imputationVariables <- sample(imputationVariables)
    }

    imp.dataset[, imputationVariables] <- missForest::prodNA(as.data.frame(xTrue[, imputationVariables]), prodNAs)
    cl <- parallel::makeCluster(clustersNumber)
    doParallel::registerDoParallel(cl)
    if (prodNAs != 0) {
      rfImp.res <- missForest::missForest(xmis = as.matrix(imp.dataset),
                                          maxiter = 50, ntree = 1000, parallelize = "forests",
                                          verbose = F, variablewise = T, decreasing = T,
                                          xtrue = as.matrix(xTrue))

      # R2 calculation
      r2.var <- numeric()
      for(impVariable in imputationVariables){
        matObs <- as.matrix(xTrue[, impVariable])
        matNA <- as.matrix(imp.dataset[, impVariable])

        observed <- matObs[which(!is.na(matObs) & is.na(matNA))]

        predicted <- as.matrix(rfImp.res$ximp[, impVariable])
        predicted <- predicted[which(!is.na(matObs) & is.na(matNA))]

        r2.var[impVariable] <- cor(observed, predicted)^2
      }

      predictivePerformance <- data.frame(Variable = imputationVariables,
                          N = length(imp.dataset[, 1]),
                          N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, impVariable])), 1]),
                          N_NA = length(imp.dataset[which(is.na(imp.dataset[, impVariable])), 1]),
                          NRMSE = rfImp.res$OOBerror[1:length(imputationVariables)],
                          R2 = r2.var,
                          Model = modelName)
    } else {
      rfImp.res <- missForest::missForest(xmis = as.matrix(imp.dataset),
                                          maxiter = 50, ntree = 1000, parallelize = "forests",
                                          verbose = F, variablewise = T, decreasing = T)

      predictivePerformance <- data.frame(Variable = imputationVariables,
                          N = length(imp.dataset[, 1]),
                          N_Obs = length(imp.dataset[which(!is.na(imp.dataset[, impVariable])), 1]),
                          N_NA = length(imp.dataset[which(is.na(imp.dataset[, impVariable])), 1]),
                          NRMSE = rfImp.res$OOBerror[1:length(imputationVariables)],
                          Model = modelName)
    }
    parallel::stopCluster(cl)

    row.names(predictivePerformance) <- NULL
    print(predictivePerformance)

    ## Results of the first round of imputation (gap filling) per iteration

    ximp <- as.data.frame(rfImp.res$ximp)
    ximp$taxon <- dataset$taxon
    predictivePerformance.all <- rbind(predictivePerformance.all, predictivePerformance)
    ximp.all <- rbind(ximp.all, ximp)
  }

  ### Results aggregation (for all iterations)

  imputationResults$ximp <- stats::aggregate(ximp.all[, -which(names(ximp) == "taxon")], by = list(ximp.all$taxon), FUN = mean) %>% dplyr::rename(taxon = Group.1)
  imputationResults$predictivePerformance <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"), which(names(predictivePerformance.all) == "Model"))],
                                                 by = list(predictivePerformance.all$Variable), FUN = mean) %>%
                                dplyr::rename(Variable = Group.1) %>%
                                dplyr::mutate(Model = modelName)

  if (IterationsNumber > 1) {
    imputationResults$ximp_sd <- stats::aggregate(ximp.all[, -which(names(ximp) == "taxon")], by = list(ximp.all$taxon), FUN = stats::sd) %>%
                                 dplyr::rename(taxon = Group.1)

    imputationResults$predictivePerformance_sd <- stats::aggregate(predictivePerformance.all[, -c(which(names(predictivePerformance.all) == "Variable"), which(names(predictivePerformance.all) == "Model"))],
                                                      by = list(predictivePerformance.all$Variable), FUN = stats::sd) %>% dplyr::rename(Variable = Group.1) %>%
                                     dplyr::mutate(Model = modelName)
    imputationResults$ximp_all_iterations <- ximp.all
    imputationResults$predictivePerformance_all_iterations <- predictivePerformance.all
  }

  ### IMPUTATION 2. Run models to impute variables following the previously determined order using imputed traits (except for the one being imputed) ####

  imp.dataset2 <-  xTrue
  ximp2 <- as.data.frame(imputationResults$ximp)

  for(impVariable in imputationVariables){

    for (n in 1:IterationsNumber) {
      ximp2[, impVariable] <- missForest::prodNA(as.data.frame(imp.dataset2[, impVariable]), prodNAs)

      ximp2 <- ximp2 %>% dplyr::select(-taxon)

      cl <- parallel::makeCluster(clustersNumber)
      doParallel::registerDoParallel(cl)
      if (prodNAs != 0) {
        rfImp.res2 <- missForest::missForest(xmis = as.matrix(ximp2),
                                            maxiter = 50, ntree = 1000, parallelize = "forests",
                                            verbose = F, variablewise = T, decreasing = T,
                                            xtrue = as.matrix(xTrue))


        # R2 calculation
        matObs <- as.matrix(xTrue[, impVariable])
        matNA <- as.matrix(ximp2[, impVariable])

        observed <- matObs[which(!is.na(matObs) & is.na(matNA))]

        predicted <- as.matrix(rfImp.res$ximp[, impVariable])
        predicted <- predicted[which(!is.na(matObs) & is.na(matNA))]

        r2.var <- cor(observed, predicted)^2

        predictivePerformance <- data.frame(Variable = impVariable,
                            N = length(ximp2[, 1]),
                            N_Obs = length(ximp2[which(!is.na(ximp2[, impVariable])), 1]),
                            N_NA = length(ximp2[which(is.na(ximp2[, impVariable])), 1]),
                            NRMSE = rfImp.res2$OOBerror[1:length(impVariable)],
                            R2 = r2.var,
                            Model = modelName)
      } else {
        rfImp.res2 <- missForest::missForest(xmis = as.matrix(ximp2),
                                             N = length(ximp2[, 1]),
                                             N_Obs = length(ximp2[which(!is.na(ximp2[, impVariable])), 1]),
                                             N_NA = length(ximp2[which(is.na(ximp2[, impVariable])), 1]),
                                            maxiter = 50, ntree = 1000, parallelize = "forests",
                                            verbose = F, variablewise = T, decreasing = T)

        predictivePerformance <- data.frame(Variable = impVariable,
                            NRMSE = rfImp.res2$OOBerror[1:length(impVariable)],
                            Model = modelName)
      }
      parallel::stopCluster(cl)

      ## Results of the second round of imputation (gap filling) per iteration

      row.names(predictivePerformance) <- NULL
      print(predictivePerformance)
      ximp2[, impVariable] <- as.data.frame(rfImp.res2$ximp[, impVariable])
      ximp2$taxon <- dataset$taxon
      predictivePerformance.all2 <- rbind(predictivePerformance.all2, predictivePerformance)
      ximp.all2 <- rbind(ximp.all2, ximp2)

      imp.dataset2[, impVariable] <- as.data.frame(rfImp.res2$ximp[, impVariable])
    }
  }

  ### Results aggregation (for all iterations)

  imputationResults$ximp2 <- stats::aggregate(ximp.all2[, -which(names(ximp.all2) == "taxon")], by = list(ximp.all2$taxon), FUN = mean) %>% dplyr::rename(taxon = Group.1)
  imputationResults$predictivePerformance2 <- stats::aggregate(predictivePerformance.all2[, -c(which(names(predictivePerformance.all2) == "Variable"), which(names(predictivePerformance.all2) == "Model"))],
                                                 by = list(predictivePerformance.all2$Variable), FUN = mean) %>%
    dplyr::rename(Variable = Group.1) %>%
    dplyr::mutate(Model = modelName)

  if (IterationsNumber > 1) {
    imputationResults$ximp_sd <- stats::aggregate(ximp.all2[, -which(names(ximp) == "taxon")], by = list(ximp.all2$taxon), FUN = stats::sd) %>%
      dplyr::rename(taxon = Group.1)

    imputationResults$predictivePerformance_sd2 <- stats::aggregate(predictivePerformance.all2[, -c(which(names(predictivePerformance.all2) == "Variable"), which(names(predictivePerformance.all2) == "Model"))],
                                                      by = list(predictivePerformance.all2$Variable), FUN = stats::sd) %>% dplyr::rename(Variable = Group.1) %>%
      dplyr::mutate(Model = modelName)
    imputationResults$ximp_all_iterations2 <- ximp.all2
    imputationResults$predictivePerformance_all_iterations2 <- predictivePerformance.all2
  }

  return(imputationResults)
}
