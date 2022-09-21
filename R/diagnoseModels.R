#' Diangose mcmcglmm model autocorrelation, convergence and effective sample size.
#'
#' @param model Mcmcglmm model
#'
#' @return Model diagnostics summary
#' @export
#'
#' @examples
diagnoseModels <- function(model, verbose = T){

  model$autocFix <- coda::autocorr.diag(model$Sol)
  model$autocRan <- coda::autocorr.diag(model$VCV)
  model$heidelFix <- coda::heidel.diag(model$Sol)
  model$heidelRan <- coda::heidel.diag(model$VCV)
  model$effSizeFix <- coda::effectiveSize(model$VCV)
  model$effSizeRan <- coda::effectiveSize(model$Sol)
  name <- model$name

  mdlDiagn <- data.frame("model" = name, "AutcorrFix" = "T", "AutcorrRan" = "T", "HeidFix" = "T", "heidRan" = "T", "effSizeFix" = "T", "effSizeRan" = "T")

  # Autocorrelation
  if(any(abs(model$autocFix[-1, ]) > 0.1)){
    if(verbose){
      cat(name, " autocFix Failed. Increase thinning. \n")
    }
    mdlDiagn$AutcorrFix <- "F"
  }
  if(any(abs(model$autocRan[-1, ]) > 0.1)){
    if(verbose){
      cat(name, " autocRan Failed. Increase thinning. \n")
    }
    cat(name, " autocRan Failed. Increase thinning. \n")
    mdlDiagn$AutcorrRan <- "F"
  }
  # Convergence (stationary)
  i <- length(summary(model)$solutions)/5
  heidelFix <- model$heidelFix[1:i, 1:3]
  if(any(heidelFix == 0)){
    if(verbose){
      cat(name, " heidelFix Failed. Increase iterations/burning. \n")
    }
    mdlDiagn$HeidFix <- "F"
  }
  i <- (length(summary(model)$Gcovariances)/4) + 1
  heidelRan <- model$heidelRan[1:i, 1:3]
  if(any(heidelRan == 0)){
    if(verbose){
      cat(name, " heidelFix Failed. Increase iterations/burning. \n")
    }
    mdlDiagn$heidRan <- "F"
  }
  # Effect size
  if(any(model$effSizeFix < 1000)){
    if(verbose){
      cat(name, " EffSizeFix Failed. Increase iterations. \n")
    }
    mdlDiagn$effSizeFix <- "F"
  }
  if(any(model$effSizeRan < 1000)){
    if(verbose){
      cat(name, " EffSizeRan Failed. Increase iterations.  \n")
    }
    mdlDiagn$effSizeRan <- "F"
  }
  return(mdlDiagn)
}
