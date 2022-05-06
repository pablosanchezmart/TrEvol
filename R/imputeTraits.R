#' Impute traits using the evolutionary (or genetic) correlations.
#'
#' @param imputationVariables (character) Names of the variables with NA where imputations will be implemented. If more than one, covariation among vimputation variables is considered.
#' @param PREDICTORS (character) Names of the variables without NAs used as predictors.
#' @param DATASET (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param PHYLOGENY (phylo) Phylogeny with tip labels contained in dataset$animal
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
imputeTraits <- function(DATASET, PHYLOGENY, correlationsTraitsResults, imputationVariables,
                             PREDICTORS = c(paste0("Phylo_axis_", 1:10)), prodNAs = 0, IterationsNumber = 10, clustersNumber = 2, FORCERUN = T){

  #### PHYLOGENETIC PRNCIPAL COMPONENTS ------------------------------------------ ####

  # Genus level phylogeny (run only once)

  if(!file.exists(paste0(outputs.dir, "phylo_eigenvectors.csv")) | isTRUE(FORCERUN)){

    mat <- ape::cophenetic.phylo(PHYLOGENY)
    mat <- as.data.frame(mat)
    pCordA <- ape::pcoa(mat)

    phyloEigenV <- as.data.frame(pCordA$vectors)
    phyloEigenV <- phyloEigenV[, c(1:50)]
    phyloEigenV$taxon <- rownames(phyloEigenV)

    imp.dataset <- merge(DATASET, phyloEigenV, by = "taxon", all.x = T)

    names(imp.df) <- stringr::str_replace_all(names(imp.df), "Axis.", "Phylo_axis_")

    utils::write.csv(imp.df, paste0(outputs.dir, "phylo_eigenvectors.csv"), row.names = F)
    print(paste0(outputs.dir, "phylo_eigenvectors.csv"))
  } else {
    # Phylogenetic principal coordinates merge
    imp.df <-  utils::read.csv(paste0(outputs.dir, "phylo_eigenvectors.csv"), header = T)
    print("loading previously calculated phylogenetic eigenvectors")
  }

  #### IMPUTATION ALGORITHM ------------------------------------------------------ ####

  phylo_order <- correlationsTraitsResults$correlation.results %>% dplyr::arrange(dplyr::desc(abs(Phylogenetic_cor))) %>% dplyr::select(Variable1, Variable2, Phylogenetic_cor)

  imputation.variables <- character()

  for(i in 1:length(phylo_order[, 1])){
    impVars <- phylo_order[i, ] %>% dplyr::select(Variable1, Variable2)
    impVars <- c(impVars$Variable1, impVars$Variable2)
    impVars <-impVars[!impVars %in% imputation.variables]
    imputation.variables <- c(imputation.variables, impVars)
  }

  imputationVariables <- imputation.variables[imputation.variables %in% imputationVariables]

  imp.df <- imp.df %>% dplyr::select(imputationVariables, PREDICTORS)
  xTrue <- imp.df

  OOBerror.all <- data.frame()
  ximp.all <- data.frame()

  for(n in 1:IterationsNumber){
    imp.df[, imputationVariables] <- missForest::prodNA(as.data.frame(imp.df[, imputationVariables]), prodNAs)

    #### PREDICTIVE MODEL ---------------------------------------------------------- #####

    cl <- parallel::makeCluster(clustersNumber)
    doParallel::registerDoParallel(cl)
    if(prodNAs != 0){
      a <- as.matrix(imp.df)
      rfImp.res <- missForest::missForest(xmis = as.matrix(imp.df), maxiter = 50, ntree = 1000, parallelize = "forests",
                                          verbose = F, variablewise = T, decreasing = T, xtrue = as.matrix(xTrue))
    }  else{
      rfImp.res <- missForest::missForest(xmis = as.matrix(imp.df), maxiter = 50, ntree = 1000, parallelize = "forests",
                                          verbose = F, variablewise = T, decreasing = T)
    }
    parallel::stopCluster(cl)
    OOBerror <- data.frame("Variable" = imputationVariables,
                           "RMSE" = sqrt(rfImp.res$OOBerror[1:length(imputationVariables)]),
                           "Mean" = apply(as.data.frame(xTrue[, imputationVariables]), MARGIN = 2, FUN = mean, na.rm = T),
                           "NRMSE" = sqrt(rfImp.res$OOBerror[1:length(imputationVariables)])/abs(apply(as.data.frame(xTrue[, imputationVariables]), MARGIN = 2, FUN = mean, na.rm = T)),
                           "Sd" = apply(as.data.frame(xTrue[, imputationVariables]), MARGIN = 2, FUN = stats::sd, na.rm = T),
                           "Model" = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))
    row.names(OOBerror) <- NULL
    print(OOBerror)

    ximp <- as.data.frame(rfImp.res$ximp)
    ximp$taxa <- DATASET$animal
    OOBerror.all <- rbind(OOBerror.all, OOBerror)
    ximp.all <- rbind(ximp.all, ximp)
  }
  rfImp <- list()

  rfImp$ximp <- stats::aggregate(ximp.all[, -which(names(ximp) == "taxa")], by = list(ximp.all$taxa), FUN = mean) %>%
    dplyr::rename(taxa = Group.1)

  rfImp$OOBerror <- stats::aggregate(OOBerror.all[, -c( which(names(OOBerror.all) == "Variable"), which(names(OOBerror.all) == "Model") )], by = list(OOBerror.all$Variable), FUN = mean) %>%
    dplyr::rename(Variable = Group.1) %>%
    dplyr::mutate("Model" = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))

  if(IterationsNumber > 1){
    rfImp$ximp_sd <- stats::aggregate(ximp.all[, -which(names(ximp) == "taxa")], by = list(ximp.all$taxa), FUN = stats::sd)%>%
      dplyr::rename(taxa = Group.1)
    rfImp$OOBerror_sd <- stats::aggregate(OOBerror.all[, -c( which(names(OOBerror.all) == "Variable"), which(names(OOBerror.all) == "Model") )], by = list(OOBerror.all$Variable), FUN = stats::sd) %>%
      dplyr::rename(Variable = Group.1) %>%
      dplyr::mutate("Model" = paste(paste(imputationVariables, collapse = ", "), paste(PREDICTORS, collapse = ", ")))

    rfImp$ximp_all_iterations <- ximp.all
    rfImp$OOBerror_all_iterations <- OOBerror.all
  }

  return(rfImp)
}
