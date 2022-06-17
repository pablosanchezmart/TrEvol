#' Impute traits using the evolutionary (or genetic) correlations.
#'
#' @param imputationVariables (character) Names of the variables with NA where imputations will be implemented. If more than one, covariation among vimputation variables is considered.
#' @param PREDICTORS (character) Names of the variables without NAs used as predictors.
#' @param DATASET (data frame) Dataset containing the variable of interest and a column named taxon describing terminal taxa of phylogeny.
#' @param PHYLOGENY (phylo) Phylogeny with tip labels contained in dataset$taxon
#' @param correlationsTraitsResults (list) Results from correlationsTraits() function.
#' @param prodNAs (numeric) Proportion of artificial NAs to be introduced in a given dataset. For evaluation puposes.
#' @param IterationsNumber (integer) Number of iterations of the imputation processs. Results are summarized.
#' @param clustersNumber (integer) Number of clusters to use in parallelization.
#' @param FORCERUN (logical) If false, previously calculated phylogenetic eigenvectors are used.
#'
#' @return
#' @export
#'
#' @examples
imputeTraits <- function(DATASET, PHYLOGENY, correlationsTraitsResults, varianceResults = NULL, imputationVariables, orderCriterium = NULL,
                             PREDICTORS = c(paste0("Phylo_axis_", 1:10)), prodNAs = 0, IterationsNumber = 10, clustersNumber = 2, FORCERUN = T){

  ### Objects to gather results ####

  OOBerror.all <- data.frame()
  ximp.all <- data.frame()

  OOBerror.all2 <- data.frame()
  ximp.all2 <- data.frame()

  imputationResults <- list()

  #### PHYLOGENETIC PRNCIPAL COMPONENTS ------------------------------------------ ####

  ### Phylogenetic principal coordinates ####

  if (!file.exists(paste0(outputs.dir, "phylo_eigenvectors.csv")) |
      isTRUE(FORCERUN)) {
    mat <- ape::cophenetic.phylo(PHYLOGENY)
    mat <- as.data.frame(mat)
    pCordA <- ape::pcoa(mat)
    phyloEigenV <- as.data.frame(pCordA$vectors)
    phyloEigenV <- phyloEigenV[DATASET$taxon, c(1:50)]
    phyloEigenV$taxon <- rownames(phyloEigenV)
    imp.dataset <- merge(DATASET, phyloEigenV, by = "taxon", all.x = T)
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

  phylo_order <- correlationsTraitsResults %>%
    dplyr::arrange(dplyr::desc(abs(Pure_coordinated_phylogenetic_conservatism))) %>%
    dplyr::select(c(Trait_1, Trait_2, orderCriterium))

  imputation.variables <- character()
  for (i in 1:length(phylo_order[, 1])) {
    impVars <- phylo_order[i, ] %>% dplyr::select(Trait_1, Trait_2)
    impVars <- c(impVars$Trait_1, impVars$Trait_2)
    impVars <- impVars[!impVars %in% imputation.variables]
    imputation.variables <- c(imputation.variables, impVars)
  }

  # Order the first two traits according to their phylogenetic variance

  firstTwoOrder <-  varianceResults %>%
    dplyr::filter(Trait %in% imputation.variables[1:2]) %>%
    dplyr::arrange(dplyr::desc(abs(Total_phylogenetic_conservatism))) %>%
    dplyr::pull(Trait)

  imputation.variables[1:2] <- firstTwoOrder

  imputationVariables <- imputation.variables[imputation.variables %in% imputationVariables]
  imp.dataset <- imp.dataset %>% dplyr::select(c(imputationVariables, PREDICTORS))
  xTrue <- imp.dataset

  ### IMPUTATION 1. Run models to impute variables following the previously determined order ####

  for (n in 1:IterationsNumber) {
    imp.dataset[, imputationVariables] <- missForest::prodNA(as.data.frame(xTrue[, imputationVariables]), prodNAs)
    cl <- parallel::makeCluster(clustersNumber)
    doParallel::registerDoParallel(cl)
    if (prodNAs != 0) {
      rfImp.res <- missForest::missForest(xmis = as.matrix(imp.dataset),
                                          maxiter = 50, ntree = 1000, parallelize = "forests",
                                          verbose = F, variablewise = T, decreasing = T,
                                          xtrue = as.matrix(xTrue))
    } else {
      rfImp.res <- missForest::missForest(xmis = as.matrix(imp.dataset),
                                          maxiter = 50, ntree = 1000, parallelize = "forests",
                                          verbose = F, variablewise = T, decreasing = T)
    }
    parallel::stopCluster(cl)

    ## Results of the first round of imputation (gap filling) per iteration

    OOBerror <- data.frame(Variable = imputationVariables,
                           RMSE = sqrt(rfImp.res$OOBerror[1:length(imputationVariables)]),
                           Mean = apply(as.data.frame(xTrue[, imputationVariables]),
                                        MARGIN = 2, FUN = mean, na.rm = T),
                           # NRMSE = sqrt(rfImp.res$OOBerror[1:length(imputationVariables)]) / abs(apply(as.data.frame(xTrue[, imputationVariables]), MARGIN = 2, FUN = mean, na.rm = T)),
                           Sd = apply(as.data.frame(xTrue[, imputationVariables]), MARGIN = 2, FUN = stats::sd, na.rm = T),
                           Model = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))
    row.names(OOBerror) <- NULL
    print(OOBerror)
    ximp <- as.data.frame(rfImp.res$ximp)
    ximp$taxon <- DATASET$taxon
    OOBerror.all <- rbind(OOBerror.all, OOBerror)
    ximp.all <- rbind(ximp.all, ximp)
  }

  ### Results aggregation (for all iterations)

  imputationResults$ximp <- stats::aggregate(ximp.all[, -which(names(ximp) == "taxon")], by = list(ximp.all$taxon), FUN = mean) %>% dplyr::rename(taxon = Group.1)
  imputationResults$OOBerror <- stats::aggregate(OOBerror.all[, -c(which(names(OOBerror.all) == "Variable"), which(names(OOBerror.all) == "Model"))],
                                                 by = list(OOBerror.all$Variable), FUN = mean) %>%
                                dplyr::rename(Variable = Group.1) %>%
                                dplyr::mutate(Model = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))

  if (IterationsNumber > 1) {
    imputationResults$ximp_sd <- stats::aggregate(ximp.all[, -which(names(ximp) == "taxon")], by = list(ximp.all$taxon), FUN = stats::sd) %>%
                                 dplyr::rename(taxon = Group.1)

    imputationResults$OOBerror_sd <- stats::aggregate(OOBerror.all[, -c(which(names(OOBerror.all) == "Variable"), which(names(OOBerror.all) == "Model"))],
                                                      by = list(OOBerror.all$Variable), FUN = stats::sd) %>% dplyr::rename(Variable = Group.1) %>%
                                     dplyr::mutate(Model = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))
    imputationResults$ximp_all_iterations <- ximp.all
    imputationResults$OOBerror_all_iterations <- OOBerror.all
  }

  ### IMPUTATION 2. Run models to impute variables following the previously determined order using imputed traits (except for the one being imputed) ####

  imp.dataset2 <-  xTrue
  ximp2 <- as.data.frame(imputationResults$ximp)

  for(impVariable in imputationVariables){

    # ximp2[, impVariable] <- imp.dataset2[, impVariable]

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
      } else {
        rfImp.res2 <- missForest::missForest(xmis = as.matrix(ximp2),
                                            maxiter = 50, ntree = 1000, parallelize = "forests",
                                            verbose = F, variablewise = T, decreasing = T)
      }
      parallel::stopCluster(cl)

      ## Results of the second round of imputation (gap filling) per iteration

      OOBerror <- data.frame(Variable = impVariable,
                             RMSE = sqrt(rfImp.res2$OOBerror[1:length(impVariable)]),
                             Mean = apply(as.data.frame(xTrue[, impVariable]),
                                          MARGIN = 2, FUN = mean, na.rm = T),
                             # NRMSE = sqrt(rfImp.res2$OOBerror[1:length(impVariable)]) / abs(apply(as.data.frame(xTrue[, impVariable]), MARGIN = 2, FUN = mean, na.rm = T)),
                             Sd = apply(as.data.frame(xTrue[, impVariable]), MARGIN = 2, FUN = stats::sd, na.rm = T),
                             Model = paste(paste(impVariable, collapse = ", "), paste(PREDICTORS, collapse = ", ")))
      row.names(OOBerror) <- NULL
      print(OOBerror)
      ximp2[, impVariable] <- as.data.frame(rfImp.res2$ximp[, impVariable])
      ximp2$taxon <- DATASET$taxon
      OOBerror.all2 <- rbind(OOBerror.all2, OOBerror)
      ximp.all2 <- rbind(ximp.all2, ximp2)

      imp.dataset2[, impVariable] <- as.data.frame(rfImp.res2$ximp[, impVariable])
    }
  }

  ### Results aggregation (for all iterations)

  imputationResults$ximp2 <- stats::aggregate(ximp.all2[, -which(names(ximp.all2) == "taxon")], by = list(ximp.all2$taxon), FUN = mean) %>% dplyr::rename(taxon = Group.1)
  imputationResults$OOBerror2 <- stats::aggregate(OOBerror.all2[, -c(which(names(OOBerror.all2) == "Variable"), which(names(OOBerror.all2) == "Model"))],
                                                 by = list(OOBerror.all2$Variable), FUN = mean) %>%
    dplyr::rename(Variable = Group.1) %>%
    dplyr::mutate(Model = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))

  if (IterationsNumber > 1) {
    imputationResults$ximp_sd <- stats::aggregate(ximp.all2[, -which(names(ximp) == "taxon")], by = list(ximp.all2$taxon), FUN = stats::sd) %>%
      dplyr::rename(taxon = Group.1)

    imputationResults$OOBerror_sd2 <- stats::aggregate(OOBerror.all2[, -c(which(names(OOBerror.all2) == "Variable"), which(names(OOBerror.all2) == "Model"))],
                                                      by = list(OOBerror.all2$Variable), FUN = stats::sd) %>% dplyr::rename(Variable = Group.1) %>%
      dplyr::mutate(Model = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))
    imputationResults$ximp_all_iterations2 <- ximp.all2
    imputationResults$OOBerror_all_iterations2 <- OOBerror.all2
  }

  return(imputationResults)
}
