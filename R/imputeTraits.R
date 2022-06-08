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
imputeTraits <- function(DATASET, PHYLOGENY, correlationsTraitsResults, imputationVariables, orderCriterium = "orderCriterium",
                             PREDICTORS = c(paste0("Phylo_axis_", 1:10)), prodNAs = 0, IterationsNumber = 10, clustersNumber = 2, FORCERUN = T){

  #### PHYLOGENETIC PRNCIPAL COMPONENTS ------------------------------------------ ####

  # Genus level phylogeny (run only once)

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
  phylo_order <- correlationsTraitsResults %>%
    dplyr::arrange(dplyr::desc(abs(orderCriterium))) %>%
    dplyr::select(Trait_1, Trait_2, orderCriterium)

  imputation.variables <- character()
  for (i in 1:length(phylo_order[, 1])) {
    impVars <- phylo_order[i, ] %>% dplyr::select(Trait_1, Trait_2)
    impVars <- c(impVars$Trait_1, impVars$Trait_2)
    impVars <- impVars[!impVars %in% imputation.variables]
    imputation.variables <- c(imputation.variables, impVars)
  }

  imputationVariables <- imputation.variables[imputation.variables %in% imputationVariables]
  imp.dataset <- imp.dataset %>% dplyr::select(imputationVariables, PREDICTORS)
  xTrue <- imp.dataset
  OOBerror.all <- data.frame()
  ximp.all <- data.frame()

  for (n in 1:IterationsNumber) {
    imp.dataset[, imputationVariables] <- missForest::prodNA(as.data.frame(imp.dataset[, imputationVariables]), prodNAs)
    cl <- parallel::makeCluster(clustersNumber)
    doParallel::registerDoParallel(cl)
    if (prodNAs != 0) {
      a <- as.matrix(imp.dataset)
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

  imputationResults <- list()
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

  return(imputationResults)
}
