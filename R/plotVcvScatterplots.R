#' Plot scatterplots showing the variance-covariancae related to phylogeny-neutral genetic structure and convergence.
#'
#' @param variable1 (character) Name of the first variable. It must be contained in dataset.
#' @param variable2 (character) Name of the second variable. It must be contained in dataset.
#' @param environmental.variables  (character) Names of the environmental.variables. They must be contained in dataset.
#' @param phylogeny (phylo) phylogeny with tip labels contained in dataset$animal
#' @param dataset (data frame) dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param model.specifications (list) Mcmcglmm models specifications. See defineModelsSpecification.
#'
#' @return
#' @export
#'
#' @examples
plotVcvScatterplots <- function(variable1, variable2, environmental.variables = NULL, dataset,
                                model.specifications, phylogeny){

  #### TOTAL RESULTS ----------------------------------------------------------- ####
  ### Data ####

  fix.frml <- paste0("cbind(", variable1, ", ", variable2, ") ~ trait-1")

  if(!is.null(environmental.variables)){
    # Complete data
    modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(variable1, variable2, environmental.variables))
    modellingData$dta <- modellingData$dta %>% dplyr::select(animal, variable1, variable2, environmental.variables)
  } else{
    modellingData <- completePhyloData(phylogeny = phylogeny, dataset = dataset, traits = c(variable1, variable2))
    modellingData$dta <- modellingData$dta %>% dplyr::select(animal, variable1, variable2)
  }

  ### Model ####

  # formula

  fix.frml <- paste0("cbind(", variable1, ", ", variable2, ") ~ trait-1")

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

  PS.v1 <- computeVariancePartition(traits = variable1, dataset = dataset, phylogeny = phylogeny, model.specifications = model.specifications, force.run = F, save = F)
  PS.v1 <- PS.v1$varianceResults$Total_phylogenetic_conservatism
  PS.v2 <- computeVariancePartition(traits = variable2, dataset = dataset, phylogeny = phylogeny, model.specifications = model.specifications, save = F)
  PS.v2 <- PS.v2$varianceResults$Total_phylogenetic_conservatism

  VCV.res <- computeCovariancePartition(traits = c(variable1, variable2), dataset = dataset, phylogeny = phylogeny, model.specifications = model.specifications, save = F)

  ## Correlations

  totalCorrlation  <-  VCV.res$covarianceResults$Total_coordination

  relativePhylogeneticCorrelation  <-  VCV.res$covarianceResults$Total_coordinated_phylogenetic_conservatism

  relativeNonPhylogeneticCorrelation <- VCV.res$covarianceResults$Total_coordinated_radiation

  # residuals

  fit <- MCMCglmm::predict.MCMCglmm(mdl, marginal = NULL)

  fit.df <- data.frame("fit_1" = fit[1:length(modellingData$dta$animal), ],
                       "fit_2" = fit[(length(modellingData$dta$animal)+1):length(fit), ])
  colnames(fit.df) <- c(variable1, variable2)

  modellingData$dta[, paste0("residuals_", variable1, ".phylogeny")] <- modellingData$dta[, variable1] - fit.df[, variable1]
  modellingData$dta[, paste0("residuals_", variable2, ".phylogeny")] <- modellingData$dta[, variable2] - fit.df[, variable2]

  modellingData$dta[, paste0("prediction_", variable1, ".phylogeny")] <- modellingData$dta[, variable1] - modellingData$dta[, paste0("residuals_", variable1, ".phylogeny")]
  modellingData$dta[, paste0("prediction_", variable2, ".phylogeny")] <- modellingData$dta[, variable2] - modellingData$dta[, paste0("residuals_", variable2, ".phylogeny")]

  ### PARTIAL RESULTS ####


  if(!is.na(environmental.variables)){

    ### Model ####

    fix.frml <- paste0("cbind(", variable1, ", ", variable2, ") ~ trait-1")
    for(predictor in environmental.variables){
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

    PS.v1.pred <- computeVariancePartition(traits = variable1, environmental.variables = environmental.variables, dataset = dataset, phylogeny = phylogeny,
                                                   model.specifications = model.specifications, force.run = T, save = F)
    PS.v1.pred <- PS.v1.pred$varianceResults$Pure_phylogenetic_conservatism

    PS.v2.pred <- computeVariancePartition(traits = variable2, environmental.variables = environmental.variables, dataset = dataset, phylogeny = phylogeny,
                                           model.specifications = model.specifications, force.run = T, save = F)
    PS.v2.pred <- PS.v2.pred$varianceResults$Pure_phylogenetic_conservatism

    VCV.pred.res <- computeCovariancePartition(traits = c(variable1, variable2), environmental.variables = environmental.variables, dataset = dataset, phylogeny = phylogeny,
                                               model.specifications = model.specifications, force.run = T, save = F)

    ## Correlations

    relativeNonPhylogeneticCorrelation.pred  <- VCV.pred.res$covarianceResults$Residual_coordination

    relativePhylogeneticCorrelation.pred  <-  VCV.pred.res$covarianceResults$Pure_coordinated_phylogenetic_conservatism

    # residuals without random effect

    fit <- MCMCglmm::predict.MCMCglmm(mdl,
                            marginal = mdl$Random$formula)

    fit.df <- data.frame("fit_1" = fit[1:length(modellingData$dta$animal), ],
                         "fit_2" = fit[(length(modellingData$dta$animal)+1):length(fit), ])
    colnames(fit.df) <- c(variable1, variable2)

    modellingData$dta[, paste0("residuals_", variable1, ".", environmental.variables)] <- modellingData$dta[, variable1] - fit.df[, variable1]
    modellingData$dta[, paste0("residuals_", variable2, ".", environmental.variables)] <- modellingData$dta[, variable2] - fit.df[, variable2]

    modellingData$dta[, paste0("prediction_", variable1, ".", environmental.variables)] <- modellingData$dta[, variable1] - modellingData$dta[, paste0("residuals_", variable1, ".", environmental.variables)]
    modellingData$dta[, paste0("prediction_", variable2, ".", environmental.variables)] <- modellingData$dta[, variable2] - modellingData$dta[, paste0("residuals_", variable2, ".", environmental.variables)]

    # residuals with random effect

    fit <- MCMCglmm::predict.MCMCglmm(mdl, marginal = NULL)

    fit.df <- data.frame("fit_1" = fit[1:length(modellingData$dta$animal), ],
                         "fit_2" = fit[(length(modellingData$dta$animal)+1):length(fit), ])
    colnames(fit.df) <- c(variable1, variable2)

    modellingData$dta[, paste0("residuals_", variable1, ".", environmental.variables, "_phylogeny")] <- modellingData$dta[, variable1] - fit.df[, variable1]
    modellingData$dta[, paste0("residuals_", variable2, ".", environmental.variables, "_phylogeny")] <- modellingData$dta[, variable2] - fit.df[, variable2]

    modellingData$dta[, paste0("prediction_", variable1, ".", environmental.variables, "_phylogeny")] <- modellingData$dta[, variable1]  -  modellingData$dta[, paste0("residuals_", variable1, ".", environmental.variables, "_phylogeny")]
    modellingData$dta[, paste0("prediction_", variable2, ".", environmental.variables, "_phylogeny")] <- modellingData$dta[, variable2] -  modellingData$dta[, paste0("residuals_", variable2, ".", environmental.variables, "_phylogeny")]

    # phylogeny - fixed effect

    modellingData$dta[, paste0("prediction_", variable1, "_phylogeny", "-environmental.variables")] <- modellingData$dta[, paste0("prediction_", variable1, ".phylogeny")]  - modellingData$dta[, paste0("prediction_", variable1, ".", environmental.variables)]
    modellingData$dta[, paste0("prediction_", variable2, "_phylogeny", "-environmental.variables")] <- modellingData$dta[, paste0("prediction_", variable2, ".phylogeny")] - modellingData$dta[, paste0("prediction_", variable2, ".", environmental.variables)]
  }

  ### PICs ####

  # force dicothomy
  modellingData$phylo <- ape::multi2di(modellingData$phylo, random=TRUE)

  modellingData$pics <- data.frame(ape::pic(modellingData$dta[, variable1], modellingData$phylo),
                                   ape::pic(modellingData$dta[, variable2], modellingData$phylo))
  names(modellingData$pics) <- c(variable1, variable2)

  rslts <- list()
  rslts$data <- modellingData

  ### TOTAL PLOTS--------------------------------------------------------------- ####

  lowX <- min(modellingData$dta[, variable2]) - abs(1.5 * max(modellingData$dta[, variable2]))
  highX <- max(modellingData$dta[, variable2]) + abs(1.5 * max(modellingData$dta[, variable2]))

  rangeX <- highX - lowX

  lowY <- min(modellingData$dta[, variable1]) - abs(1.5 * max(modellingData$dta[, variable1]))
  highY <- max(modellingData$dta[, variable1])  + abs(1.5 * max(modellingData$dta[, variable1]))

  rangeY <- highY - lowY

  variable1_axis <- stringr::str_replace_all(variable1, "_", " ")
  variable1_axis <- stringr::str_remove_all(variable1_axis, "ln")

  variable2_axis <- stringr::str_replace_all(variable2, "_", " ")
  variable2_axis <- stringr::str_remove_all(variable2_axis, "ln")

  ### Data correlation ####

  pCor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, variable1], x = modellingData$dta[, variable2])) + ggplot2::geom_point(alpha = 0.5)
  pCor <- pCor + ggplot2::stat_ellipse(level = 0.95)
  pCor <- pCor +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(mean(totalCorrlation), 2)),
                           x = Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")

  pCor <- pCor + ggplot2::xlab(variable2_axis) + ggplot2::ylab(variable1_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
  pCor <- pCor + ggplot2::ggtitle("Traits")


  rslts$pCor <- pCor


  ### PICs correlation ####

  pCor.pic <- ggplot2::ggplot(modellingData$pics, ggplot2::aes(y = modellingData$pics[, variable1], x = modellingData$pics[, variable2])) + ggplot2::geom_point(alpha = 0.5)
  pCor.pic <- pCor.pic + ggplot2::stat_ellipse(level = 0.95)
  pCor.pic <- pCor.pic +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(stats::cor(modellingData$pics[, variable1], modellingData$pics[, variable2]), 2)),
                                   x = Inf, y = -Inf, hjust = 1, vjust = 0)
  pCor.pic <- pCor.pic + ggplot2::xlab(paste0("PIC ", variable2_axis)) + ggplot2::ylab(paste0("PIC ", variable1_axis)) +
    ggplot2::xlim((mean(modellingData$pics[, variable2]) - rangeX/2), (mean(modellingData$pics[, variable2]) + rangeX/2)) +
    ggplot2::ylim((mean(modellingData$pics[, variable1]) - rangeY/2), (mean(modellingData$pics[, variable1]) + rangeY/2))
  pCor.pic <- pCor.pic + ggplot2::ggtitle("PICs")

  rslts$pCor.pic <- pCor.pic


  ### Convergence correlation (residuals . phylogeny) ####

  v1.res <- paste0("residuals_", variable1, ".phylogeny")
  v2.res <- paste0("residuals_", variable2, ".phylogeny")

  v1.res_axis <- paste0("res ", variable1_axis, ".phylo")
  v2.res_axis <- paste0("res ", variable2_axis, ".phylo")

  pCor.phylo <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res], x = modellingData$dta[, v2.res])) + ggplot2::geom_point(alpha = 0.5)
  pCor.phylo <- pCor.phylo + ggplot2::stat_ellipse(level = 0.95)
  pCor.phylo <- pCor.phylo + ggplot2::xlab(v2.res_axis) + ggplot2::ylab(v1.res_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
  pCor.phylo <- pCor.phylo +  ggplot2::annotate(geom = "text", label = paste0("relNP_r = ", round(mean(relativeNonPhylogeneticCorrelation), 2)),
                                       x = Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
  pCor.phylo <- pCor.phylo + ggplot2::ggtitle("Evolutionary convergence")
  rslts$pCor.phylo <- pCor.phylo


  ### Phylogenetic correlation (Predictions phylogeny) ####

  v1.pred <- paste0("prediction_", variable1, ".phylogeny")
  v2.pred <- paste0("prediction_", variable2, ".phylogeny")

  v1.pred_axis <- paste0("pred ", variable1_axis, ".phylo")
  v2.pred_axis <- paste0("pred ", variable2_axis, ".phylo")

  pPhyloCor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.pred], x = modellingData$dta[, v2.pred])) + ggplot2::geom_point(alpha = 0.5)
  pPhyloCor <- pPhyloCor + ggplot2::stat_ellipse(level = 0.95)
  pPhyloCor <- pPhyloCor + ggplot2::annotate(geom = "text", label = paste0("PS ", variable1_axis, " = ", round(mean(PS.v1), 2),
                                                                  "\nPS ", variable2_axis, " = ", round(mean(PS.v2), 2)),
                                    x = -Inf, y = Inf, hjust = 0, vjust = 1)
  pPhyloCor <- pPhyloCor + ggplot2::xlab(v2.pred_axis) + ggplot2::ylab(v1.pred_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
  pPhyloCor <- pPhyloCor +  ggplot2::annotate(geom = "text", label = paste0("relP_r = ", round(mean(relativePhylogeneticCorrelation), 2)),
                                     x =Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
  pPhyloCor <- pPhyloCor + ggplot2::ggtitle("Phylogenetic conservatism")

  rslts$pPhyloCor <- pPhyloCor


  ### PARTIAL PLOTS ------------------------------------------------------------ ####

  if(!is.na(environmental.variables)) {

    # relationship with predictor

    ### Correlation with predictor ####

    if(length(environmental.variables) == 1){

      predictor_axis <- stringr::str_replace_all(environmental.variables, "_", " ")

      v1.predictor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, variable1], x = modellingData$dta[, environmental.variables])) + ggplot2::geom_point(alpha = 0.5)
      v1.predictor <- v1.predictor + ggplot2::stat_ellipse(level = 0.95)
      v1.predictor <- v1.predictor +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(stats::cor(modellingData$dta[, variable1], modellingData$dta[, environmental.variables]), 2)),
                                               x = Inf, y = -Inf, hjust = 1, vjust = 0)

      v1.predictor <- v1.predictor + ggplot2::xlab(environmental.variables) + ggplot2::ylab(variable1_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
      v1.predictor <- v1.predictor + ggplot2::ggtitle(paste0(variable1_axis, " vs\n", predictor_axis))

      rslts$v1.predictor <- v1.predictor

      v2.predictor <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, variable2], x = modellingData$dta[, environmental.variables])) + ggplot2::geom_point(alpha = 0.5)
      v2.predictor <- v2.predictor + ggplot2::stat_ellipse(level = 0.95)
      v2.predictor <- v2.predictor +  ggplot2::annotate(geom = "text", label = paste0("r = ", round(stats::cor(modellingData$dta[, variable2], modellingData$dta[, environmental.variables]), 2)),
                                               x = Inf, y = -Inf, hjust = 1, vjust = 0)

      v2.predictor <- v2.predictor + ggplot2::xlab(predictor_axis) + ggplot2::ylab(variable2_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
      v2.predictor <- v2.predictor + ggplot2::ggtitle(paste0(variable2_axis, " vs\n", predictor_axis))

      rslts$v2.predictor <- v2.predictor
    }

    ### Predictor explanation power (residuals. environmental.variables) ####

    v1.res.pred <- paste0("residuals_", variable1, ".", environmental.variables)
    v2.res.pred <- paste0("residuals_", variable2, ".", environmental.variables)

    v1.res.pred_axis <- paste0("res ", variable1_axis, ".", predictor_axis)
    v2.res.pred_axis <- paste0("res ", variable2_axis, ".", predictor_axis)

    pCor.pred <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res.pred], x = modellingData$dta[, v2.res.pred])) + ggplot2::geom_point(alpha = 0.5)
    pCor.pred <- pCor.pred + ggplot2::stat_ellipse(level = 0.95)

    pCor.pred <- pCor.pred + ggplot2::xlab(v1.res.pred_axis) + ggplot2::ylab(v2.res.pred_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
    pCor.pred <- pCor.pred

    rslts$pCor.pred <- pCor.pred

    ### Partial evoltuionary convergence (residuals. environmental.variables & phylogeny) ####

    v1.res.pred.phylo <- paste0("residuals_", variable1, ".", environmental.variables, "_phylogeny")
    v2.res.pred.phylo <- paste0("residuals_", variable2, ".", environmental.variables, "_phylogeny")

    v1.res.pred.phylo_axis <- paste0("res ", variable1_axis, ".", predictor_axis, "& phylo")
    v2.res.pred.phylo_axis <- paste0("res ", variable2_axis, ".", predictor_axis, "& phylo")

    pCor.pred.phylo <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res.pred.phylo], x = modellingData$dta[, v2.res.pred.phylo])) + ggplot2::geom_point(alpha = 0.5)
    pCor.pred.phylo <- pCor.pred.phylo + ggplot2::stat_ellipse(level = 0.95)
    pCor.pred.phylo <- pCor.pred.phylo + ggplot2::xlab(v2.res.pred.phylo_axis) + ggplot2::ylab(v1.res.pred.phylo_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
    pCor.pred.phylo <- pCor.pred.phylo +  ggplot2::annotate(geom = "text", label = paste0("relNP_r = ", round(mean(relativeNonPhylogeneticCorrelation.pred), 2)),
                                                   x =Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
    pCor.pred.phylo <- pCor.pred.phylo + ggplot2::ggtitle(paste0("Partial evolutionary convergence |\n", predictor_axis))

    rslts$pCor.pred.phylo <- pCor.pred.phylo

    ### Partial phylogenetic consevatism (prediction. phylogeny - environmental.variables) ####

    v1.res.phylo_pred <- paste0("prediction_", variable1, "_phylogeny", "-environmental.variables")
    v2.resphylo_pred <- paste0("prediction_", variable2, "_phylogeny", "-environmental.variables")

    v1.resphylo_pred_axis <- paste0("pred ", variable1_axis, ".", " phylo - ", predictor_axis)
    v2.resphylo_pred_axis <- paste0("pred ", variable2_axis, ".", " phylo - ", predictor_axis)

    pPhyloCor.pred <- ggplot2::ggplot(modellingData$dta, ggplot2::aes(y = modellingData$dta[, v1.res.phylo_pred], x = modellingData$dta[, v2.resphylo_pred])) + ggplot2::geom_point(alpha = 0.5)
    pPhyloCor.pred <- pPhyloCor.pred + ggplot2::stat_ellipse(level = 0.95)
    pPhyloCor.pred <- pPhyloCor.pred + ggplot2::xlab(v2.resphylo_pred_axis) + ggplot2::ylab(v1.resphylo_pred_axis) + ggplot2::xlim(lowX, highX) + ggplot2::ylim(lowY, highY)
    pPhyloCor.pred <- pPhyloCor.pred +  ggplot2::annotate(geom = "text", label = paste0("PS ", variable1_axis, " = ", round(mean(PS.v1.pred), 2),
                                                                               "\nPS ", variable2_axis, " = ", round(mean(PS.v2.pred), 2)),
                                                 x = -Inf, y = Inf, hjust = 0, vjust = 1)
    pPhyloCor.pred <- pPhyloCor.pred +  ggplot2::annotate(geom = "text", label = paste0("relP_r = ", round(mean(relativePhylogeneticCorrelation.pred), 2)),
                                                 x =Inf, y = -Inf, hjust = 1, vjust = 0, colour = "red")
    pPhyloCor.pred <- pPhyloCor.pred + ggplot2::ggtitle(paste0("Partial phylogenetic conservatism |\n", predictor_axis))

    rslts$pPhyloCor.pred <- pPhyloCor.pred
  }
  return(rslts)
}
