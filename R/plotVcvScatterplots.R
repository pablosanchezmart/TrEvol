#' Plot scatterplots showing the variance-covariancae related to phylogeny-neutral genetic structure and convergence.
#'
#' @param VARIABLE1 (character) Name of the first variable. It must be contained in dataset.
#' @param VARIABLE2 (character) Name of the second variable. It must be contained in dataset.
#' @param PREDICTORS  (character) Names of the predictors. They must be contained in dataset.
#' @param PHYLOGENY (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param DATASET (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param MODEL.SPECIFICATIONS (list) Mcmcglmm models specifications. See defineModelsSpecification.
#'
#' @return
#' @export
#'
#' @examples
plotVcvScatterplots <- function(VARIABLE1 = "ln_Ks", VARIABLE2 = "ln_negP50", PREDICTORS, DATASET,
                                MODEL.SPECIFICATIONS, PHYLOGENY){

  #### TOTAL RESULTS ----------------------------------------------------------- ####
  ### Data ####

  fix.frml <- paste0("cbind(", VARIABLE1, ", ", VARIABLE2, ") ~ trait-1")

  if(!is.na(PREDICTORS)){
    # Complete data
    modellingData <- completePhyloData(phylogeny = PHYLOGENY, dataset = DATASET, traits = c(VARIABLE1, VARIABLE2, PREDICTORS))
    modellingData$dta <- modellingData$dta %>% dplyr::select(animal, VARIABLE1, VARIABLE2, PREDICTORS)
  } else{
    modellingData <- completePhyloData(phylogeny = PHYLOGENY, dataset = DATASET, traits = c(VARIABLE1, VARIABLE2))
    modellingData$dta <- modellingData$dta %>% dplyr::select(animal, VARIABLE1, VARIABLE2)
  }

  ### Model ####

  # formula

  fix.frml <- paste0("cbind(", VARIABLE1, ", ", VARIABLE2, ") ~ trait-1")

  # model

  mdl <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                  random = ~ us(trait):animal, rcov = ~us(trait):units,
                  data= modellingData$dta, pedigree = modellingData$phylo,
                  family = c("gaussian", "gaussian"),
                  prior = model.specifications$multiresponse_prior,
                  nitt = model.specifications$number_interations,
                  burnin = model.specifications$burning_iterations,
                  thin = model.specifications$thinning_iterations,
                  verbose = F,
                  pr = T)
  mdl$name <- fix.frml


  ### Variance-covariance calculations ####

  PS.v1 <- computePhylogeneticSignal(variable = VARIABLE1, dataset = DATASET, phylogeny = PHYLOGENY, model.specifications = MODEL.SPECIFICATIONS)
  PS.v1 <- PS.v1$phyloSignal$Wlambda
  PS.v2 <- computePhylogeneticSignal(variable = VARIABLE2, dataset = DATASET, phylogeny = PHYLOGENY, model.specifications = MODEL.SPECIFICATIONS)
  PS.v2 <- PS.v2$phyloSignal$Wlambda

  VCV.res <- computeCorrelations(variable1 = VARIABLE1, variable2 = VARIABLE2, dataset = DATASET, phylogeny = PHYLOGENY, model.specifications = MODEL.SPECIFICATIONS)

  ## Correlations

  totalCorrlation  <-  VCV.res$corrrelationsSummary$Total_cor

  phylogeneticCorrelation  <-  VCV.res$corrrelationsSummary$Phylogenetic_cor

  relativePhylogeneticCorrelation <-  VCV.res$corrrelationsSummary$Relative_phylogenetic_cor

  nonPhylogeneticCorrelation  <-  VCV.res$corrrelationsSummary$Convergent_cor

  relativeNonPhylogeneticCorrelation <- VCV.res$corrrelationsSummary$Relative_convergent_cor

  # residuals

  fit <- MCMCglmm::predict.MCMCglmm(mdl, marginal = NULL)

  fit.df <- data.frame("fit_1" = fit[1:length(modellingData$dta$animal), ],
                       "fit_2" = fit[(length(modellingData$dta$animal)+1):length(fit), ])
  colnames(fit.df) <- c(VARIABLE1, VARIABLE2)

  modellingData$dta[, paste0("residuals_", VARIABLE1, ".phylogeny")] <- modellingData$dta[, VARIABLE1] - fit.df[, VARIABLE1]
  modellingData$dta[, paste0("residuals_", VARIABLE2, ".phylogeny")] <- modellingData$dta[, VARIABLE2] - fit.df[, VARIABLE2]

  modellingData$dta[, paste0("prediction_", VARIABLE1, ".phylogeny")] <- modellingData$dta[, VARIABLE1] - modellingData$dta[, paste0("residuals_", VARIABLE1, ".phylogeny")]
  modellingData$dta[, paste0("prediction_", VARIABLE2, ".phylogeny")] <- modellingData$dta[, VARIABLE2] - modellingData$dta[, paste0("residuals_", VARIABLE2, ".phylogeny")]

  ### PARTIAL RESULTS ####


  if(!is.na(PREDICTORS)){

    ### Model ####

    fix.frml <- paste0("cbind(", VARIABLE1, ", ", VARIABLE2, ") ~ trait-1")
    for(predictor in PREDICTORS){
      fix.frml <- paste0(fix.frml, " + trait:", predictor)
    }

    mdl <- MCMCglmm::MCMCglmm(fixed = stats::as.formula(fix.frml),
                    random = ~ us(trait):animal, rcov = ~us(trait):units,
                    data= modellingData$dta, pedigree = modellingData$phylo,
                    family = c("gaussian", "gaussian"),
                    prior = model.specifications$multiresponse_prior,
                    nitt = model.specifications$number_interations,
                    burnin = model.specifications$burning_iterations,
                    thin = model.specifications$thinning_iterations,
                    verbose = F,
                    pr = T)
    mdl$name <- fix.frml

    ### Variance-covariance calculations ####

    PS.v1.pred <- computePartialPhylogeneticSignal(variable = VARIABLE1, predictors = PREDICTORS, dataset = DATASET, phylogeny = PHYLOGENY,
                                                   model.specifications = MODEL.SPECIFICATIONS)
    PS.v1.pred <- PS.v1.pred$phyloSignal$Wlambda
    PS.v2.pred <- computePartialPhylogeneticSignal(variable = VARIABLE2, predictors = PREDICTORS, dataset = DATASET, phylogeny = PHYLOGENY,
                                                   model.specifications = MODEL.SPECIFICATIONS)
    PS.v2.pred <- PS.v2.pred$phyloSignal$Wlambda

    VCV.pred.res <- computePartialCorrelations(variable1 = VARIABLE1, variable2 = VARIABLE2, predictors = PREDICTORS, dataset = DATASET, phylogeny = PHYLOGENY,
                                               model.specifications = MODEL.SPECIFICATIONS)

    ## Correlations

    totalCorrlation.pred  <- VCV.pred.res$corrrelationsSummary$Total_cor

    phylogeneticCorrelation.pred  <-  VCV.pred.res$corrrelationsSummary$Phylogenetic_cor

    relativePhylogeneticCorrelation.pred <-  VCV.pred.res$corrrelationsSummary$Relative_phylogenetic_cor

    nonPhylogeneticCorrelation.pred  <-  VCV.pred.res$corrrelationsSummary$Convergent_cor

    relativeNonPhylogeneticCorrelation.pred <-  VCV.pred.res$corrrelationsSummary$Relative_convergent_cor


    # residuals without random effect

    fit <- MCMCglmm::predict.MCMCglmm(mdl,
                            marginal = mdl$Random$formula)

    fit.df <- data.frame("fit_1" = fit[1:length(modellingData$dta$animal), ],
                         "fit_2" = fit[(length(modellingData$dta$animal)+1):length(fit), ])
    colnames(fit.df) <- c(variable1, variable2)

    modellingData$dta[, paste0("residuals_", VARIABLE1, ".", PREDICTORS)] <- modellingData$dta[, VARIABLE1] - fit.df[, VARIABLE1]
    modellingData$dta[, paste0("residuals_", VARIABLE2, ".", PREDICTORS)] <- modellingData$dta[, VARIABLE2] - fit.df[, VARIABLE2]

    modellingData$dta[, paste0("prediction_", VARIABLE1, ".", PREDICTORS)] <- modellingData$dta[, VARIABLE1] - modellingData$dta[, paste0("residuals_", VARIABLE1, ".", PREDICTORS)]
    modellingData$dta[, paste0("prediction_", VARIABLE2, ".", PREDICTORS)] <- modellingData$dta[, VARIABLE2] - modellingData$dta[, paste0("residuals_", VARIABLE2, ".", PREDICTORS)]

    # residuals with random effect

    fit <- MCMCglmm::predict.MCMCglmm(mdl,
                            marginal = NULL)

    fit.df <- data.frame("fit_1" = fit[1:length(modellingData$dta$animal), ],
                         "fit_2" = fit[(length(modellingData$dta$animal)+1):length(fit), ])
    colnames(fit.df) <- c(VARIABLE1, VARIABLE2)

    modellingData$dta[, paste0("residuals_", VARIABLE1, ".", PREDICTORS, "_phylogeny")] <- modellingData$dta[, VARIABLE1] - fit.df[, VARIABLE1]
    modellingData$dta[, paste0("residuals_", VARIABLE2, ".", PREDICTORS, "_phylogeny")] <- modellingData$dta[, VARIABLE2] - fit.df[, VARIABLE2]

    modellingData$dta[, paste0("prediction_", VARIABLE1, ".", PREDICTORS, "_phylogeny")] <- modellingData$dta[, VARIABLE1]  -  modellingData$dta[, paste0("residuals_", VARIABLE1, ".", PREDICTORS, "_phylogeny")]
    modellingData$dta[, paste0("prediction_", VARIABLE2, ".", PREDICTORS, "_phylogeny")] <- modellingData$dta[, VARIABLE2] -  modellingData$dta[, paste0("residuals_", VARIABLE2, ".", PREDICTORS, "_phylogeny")]

    # phylogeny - fixed effect

    modellingData$dta[, paste0("prediction_", VARIABLE1, "_phylogeny", "-predictors")] <- modellingData$dta[, paste0("prediction_", VARIABLE1, ".phylogeny")]  - modellingData$dta[, paste0("prediction_", VARIABLE1, ".", PREDICTORS)]
    modellingData$dta[, paste0("prediction_", VARIABLE2, "_phylogeny", "-predictors")] <- modellingData$dta[, paste0("prediction_", VARIABLE2, ".phylogeny")] - modellingData$dta[, paste0("prediction_", VARIABLE2, ".", PREDICTORS)]
  }

  ### PICs ####

  # force dicothomy
  modellingData$phylo <- ape::multi2di(modellingData$phylo, random=TRUE)

  modellingData$pics <- data.frame(ape::pic(modellingData$dta[, VARIABLE1], modellingData$phylo),
                                   ape::pic(modellingData$dta[, VARIABLE2], modellingData$phylo))
  names(modellingData$pics) <- c(VARIABLE1, VARIABLE2)

  rslts <- list()
  rslts$data <- modellingData

  ### TOTAL PLOTS--------------------------------------------------------------- ####

  lowX <- min(modellingData$dta[, VARIABLE2]) - abs(1.5 * max(modellingData$dta[, VARIABLE2]))
  highX <- max(modellingData$dta[, VARIABLE2]) + abs(1.5 * max(modellingData$dta[, VARIABLE2]))

  rangeX <- highX - lowX

  lowY <- min(modellingData$dta[, VARIABLE1]) - abs(1.5 * max(modellingData$dta[, VARIABLE1]))
  highY <- max(modellingData$dta[, VARIABLE1])  + abs(1.5 * max(modellingData$dta[, VARIABLE1]))

  rangeY <- highY - lowY

  variable1_axis <- stringr::str_replace_all(variable1, "_", " ")
  variable1_axis <- stringr::str_remove_all(variable1_axis, "ln")

  variable2_axis <- stringr::str_replace_all(variable2, "_", " ")
  variable2_axis <- stringr::str_remove_all(variable2_axis, "ln")

  ### Data correlation ####

  pCor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, VARIABLE1], x = modellingData$dta[, VARIABLE2])) + ggplot2::geom_point(alpha = 0.5)
  pCor <- pCor + ggplot2::stat_ellipse(level = 0.95)
  pCor <- pCor +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(mean(totalCorrlation), 2)),
                           x = Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")

  pCor <- pCor + ggplot2::xlab(variable2_axis) + ggplot2::ylab(variable1_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
  pCor <- pCor + ggplot2::ggtitle("Traits")


  rslts$pCor <- pCor


  ### PICs correlation ####

  pCor.pic <- ggplot2::ggplot(modellingData$pics, ggplot2::aes(y = modellingData$pics[, VARIABLE1], x = modellingData$pics[, VARIABLE2])) + ggplot2::geom_point(alpha = 0.5)
  pCor.pic <- pCor.pic + ggplot2::stat_ellipse(level = 0.95)
  pCor.pic <- pCor.pic +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(stats::cor(modellingData$pics[, VARIABLE1], modellingData$pics[, VARIABLE2]), 2)),
                                   x = Inf, y = -Inf, hjust = 1, vjust = 0)
  pCor.pic <- pCor.pic + ggplot2::xlab(paste0("PIC ", variable2_axis)) + ggplot2::ylab(paste0("PIC ", variable1_axis)) +
    ggplot2::xlim((mean(modellingData$pics[, VARIABLE2]) - rangeX/2), (mean(modellingData$pics[, VARIABLE2]) + rangeX/2)) +
    ggplot2::ylim((mean(modellingData$pics[, VARIABLE1]) - rangeY/2), (mean(modellingData$pics[, VARIABLE1]) + rangeY/2))
  pCor.pic <- pCor.pic + ggplot2::ggtitle("PICs")

  rslts$pCor.pic <- pCor.pic


  ### Convergence correlation (residuals . phylogeny) ####

  v1.res <- paste0("residuals_", VARIABLE1, ".phylogeny")
  v2.res <- paste0("residuals_", VARIABLE2, ".phylogeny")

  v1.res_axis <- paste0("res ", variable1_axis, ".phylo")
  v2.res_axis <- paste0("res ", variable2_axis, ".phylo")

  pCor.phylo <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res], x = modellingData$dta[, v2.res])) + ggplot2::geom_point(alpha = 0.5)
  pCor.phylo <- pCor.phylo + ggplot2::stat_ellipse(level = 0.95)
  pCor.phylo <- pCor.phylo + ggplot2::xlab(v2.res_axis) + ggplot2::ylab(v1.res_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
  pCor.phylo <- pCor.phylo +  ggplot2::annotate(geom = "text", label = paste0("absNP_r = ", round(mean(nonPhylogeneticCorrelation), 2)),
                                       x = -Inf, y = -Inf, hjust = 0, vjust = 0, colour = "red")
  pCor.phylo <- pCor.phylo +  ggplot2::annotate(geom = "text", label = paste0("relNP_r = ", round(mean(relativeNonPhylogeneticCorrelation), 2)),
                                       x = Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
  pCor.phylo <- pCor.phylo + ggplot2::ggtitle("Evolutionary convergence")
  rslts$pCor.phylo <- pCor.phylo


  ### Phylogenetic correlation (Predictions phylogeny) ####

  v1.pred <- paste0("prediction_", VARIABLE1, ".phylogeny")
  v2.pred <- paste0("prediction_", VARIABLE2, ".phylogeny")

  v1.pred_axis <- paste0("pred ", variable1_axis, ".phylo")
  v2.pred_axis <- paste0("pred ", variable2_axis, ".phylo")

  pPhyloCor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.pred], x = modellingData$dta[, v2.pred])) + ggplot2::geom_point(alpha = 0.5)
  pPhyloCor <- pPhyloCor + ggplot2::stat_ellipse(level = 0.95)
  pPhyloCor <- pPhyloCor + ggplot2::annotate(geom = "text", label = paste0("PS ", variable1_axis, " = ", round(mean(PS.v1), 2),
                                                                  "\nPS ", variable2_axis, " = ", round(mean(PS.v2), 2)),
                                    x = -Inf, y = Inf, hjust = 0, vjust = 1)
  pPhyloCor <- pPhyloCor + ggplot2::xlab(v2.pred_axis) + ggplot2::ylab(v1.pred_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
  pPhyloCor <- pPhyloCor +  ggplot2::annotate(geom = "text", label = paste0("absP_r = ", round(mean(phylogeneticCorrelation), 2)),
                                     x = -Inf, y = -Inf, hjust = 0, vjust = 0, colour = "red")
  pPhyloCor <- pPhyloCor +  ggplot2::annotate(geom = "text", label = paste0("relP_r = ", round(mean(relativePhylogeneticCorrelation), 2)),
                                     x =Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
  pPhyloCor <- pPhyloCor + ggplot2::ggtitle("Phylogenetic conservatism")

  rslts$pPhyloCor <- pPhyloCor


  ### PARTIAL PLOTS ------------------------------------------------------------ ####

  if(!is.na(PREDICTORS)) {

    # relationship with predictor

    ### Correlation with predictor ####

    if(length(PREDICTORS) == 1){

      predictor_axis <- stringr::str_replace_all(PREDICTORS, "_", " ")

      v1.predictor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, VARIABLE1], x = modellingData$dta[, PREDICTORS])) + ggplot2::geom_point(alpha = 0.5)
      v1.predictor <- v1.predictor + ggplot2::stat_ellipse(level = 0.95)
      v1.predictor <- v1.predictor +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(stats::cor(modellingData$dta[, VARIABLE1], modellingData$dta[, PREDICTORS]), 2)),
                                               x = Inf, y = -Inf, hjust = 1, vjust = 0)

      v1.predictor <- v1.predictor + ggplot2::xlab(PREDICTORS) + ggplot2::ylab(variable1_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
      v1.predictor <- v1.predictor

      rslts$v1.predictor <- v1.predictor

      v2.predictor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, VARIABLE2], x = modellingData$dta[, PREDICTORS])) + ggplot2::geom_point(alpha = 0.5)
      v2.predictor <- v2.predictor + ggplot2::stat_ellipse(level = 0.95)
      v2.predictor <- v2.predictor +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(stats::cor(modellingData$dta[, VARIABLE2], modellingData$dta[, PREDICTORS]), 2)),
                                               x = Inf, y = -Inf, hjust = 1, vjust = 0)

      v2.predictor <- v2.predictor + ggplot2::xlab(predictor_axis) + ggplot2::ylab(variable2_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
      v2.predictor <- v2.predictor

      rslts$v2.predictor <- v2.predictor
    }

    ### Predictor explanation power (residuals. predictors) ####

    v1.res.pred <- paste0("residuals_", VARIABLE1, ".", PREDICTORS)
    v2.res.pred <- paste0("residuals_", VARIABLE2, ".", PREDICTORS)

    v1.res.pred_axis <- paste0("res ", variable1_axis, ".", predictor_axis)
    v2.res.pred_axis <- paste0("res ", variable2_axis, ".", predictor_axis)

    pCor.pred <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res.pred], x = modellingData$dta[, v2.res.pred])) + ggplot2::geom_point(alpha = 0.5)
    pCor.pred <- pCor.pred + ggplot2::stat_ellipse(level = 0.95)

    pCor.pred <- pCor.pred + ggplot2::xlab(v1.res.pred_axis) + ggplot2::ylab(v2.res.pred_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
    pCor.pred <- pCor.pred

    rslts$pCor.pred <- pCor.pred

    ### Partial evoltuionary convergence (residuals. predictors & phylogeny) ####

    v1.res.pred.phylo <- paste0("residuals_", VARIABLE1, ".", PREDICTORS, "_phylogeny")
    v2.res.pred.phylo <- paste0("residuals_", VARIABLE2, ".", PREDICTORS, "_phylogeny")

    v1.res.pred.phylo_axis <- paste0("res ", variable1_axis, ".", predictor_axis, "& phylo")
    v2.res.pred.phylo_axis <- paste0("res ", variable2_axis, ".", predictor_axis, "& phylo")

    pCor.pred.phylo <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res.pred.phylo], x = modellingData$dta[, v2.res.pred.phylo])) + ggplot2::geom_point(alpha = 0.5)
    pCor.pred.phylo <- pCor.pred.phylo + ggplot2::stat_ellipse(level = 0.95)
    pCor.pred.phylo <- pCor.pred.phylo + ggplot2::xlab(v2.res.pred.phylo_axis) + ggplot2::ylab(v1.res.pred.phylo_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
    pCor.pred.phylo <- pCor.pred.phylo +  ggplot2::annotate(geom = "text", label = paste0("absNP_r = ", round(mean(nonPhylogeneticCorrelation.pred), 2)),
                                                   x = -Inf, y = -Inf, hjust = 0, vjust = 0, colour = "red")
    pCor.pred.phylo <- pCor.pred.phylo +  ggplot2::annotate(geom = "text", label = paste0("relNP_r = ", round(mean(relativeNonPhylogeneticCorrelation.pred), 2)),
                                                   x =Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
    pCor.pred.phylo <- pCor.pred.phylo + ggplot2::ggtitle(paste0("Partial evolutionary convergence |\n", predictor_axis))

    rslts$pCor.pred.phylo <- pCor.pred.phylo

    ### Partial phylogenetic consevatism (prediction. phylogeny - predictors) ####

    v1.res.phylo_pred <- paste0("prediction_", VARIABLE1, "_phylogeny", "-predictors")
    v2.resphylo_pred <- paste0("prediction_", VARIABLE2, "_phylogeny", "-predictors")

    v1.resphylo_pred_axis <- paste0("pred ", variable1_axis, ".", " phylo - ", predictor_axis)
    v2.resphylo_pred_axis <- paste0("pred ", variable2_axis, ".", " phylo - ", predictor_axis)

    pPhyloCor.pred <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res.phylo_pred], x = modellingData$dta[, v2.resphylo_pred])) + ggplot2::geom_point(alpha = 0.5)
    pPhyloCor.pred <- pPhyloCor.pred + ggplot2::stat_ellipse(level = 0.95)
    pPhyloCor.pred <- pPhyloCor.pred + ggplot2::xlab(v2.resphylo_pred_axis) + ggplot2::ylab(v1.resphylo_pred_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
    pPhyloCor.pred <- pPhyloCor.pred +  ggplot2::annotate(geom = "text", label = paste0("PS ", variable1_axis, " = ", round(mean(PS.v1.pred), 2),
                                                                               "\nPS ", variable2_axis, " = ", round(mean(PS.v2.pred), 2)),
                                                 x = -Inf, y = Inf, hjust = 0, vjust = 1)
    pPhyloCor.pred <- pPhyloCor.pred +  ggplot2::annotate(geom = "text", label = paste0("absP_r = ", round(mean(phylogeneticCorrelation.pred), 2)),
                                                 x = -Inf, y = -Inf, hjust = 0, vjust = 0, colour = "red")
    pPhyloCor.pred <- pPhyloCor.pred +  ggplot2::annotate(geom = "text", label = paste0("relP_r = ", round(mean(relativePhylogeneticCorrelation.pred), 2)),
                                                 x =Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
    pPhyloCor.pred <- pPhyloCor.pred + ggplot2::ggtitle(paste0("Partial phylogenetic conservatism |\n", predictor_axis))

    rslts$pPhyloCor.pred <- pPhyloCor.pred
  }
  return(rslts)
}
