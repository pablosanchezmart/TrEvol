#' Plot pairwise distances using distanace matrices
#'
#' @param matDistX (matrix) First distance matrix, distances will be displayed in the X axis.
#' @param matDistY (matrix) Second distance matrix, distances will be displayed in the Y axis.
#' @param matDistXLab (character) X axis label
#' @param matDistYLab (character) Y axis label
#'
#' @return
#' @export
#'
#' @examples
plotPairwiseDistances <- function(matDistX = NULL, matDistY = NULL, matDistXLab = "Genetic distance", matDistYLab = "Chemical distance"){

  # Packages needed

  require(graph4lg)
  require(ggplot2)

  # Function to make matrix symmetric

  makeSymm <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }


  # Make matrix symmetric

  matDistX_sym <- makeSymm(matDistX)
  matDistY_sym <- makeSymm(matDistY)

  ### Select populations with data ####
  populationsX <- rownames(matDistX_sym)
  populationsY <- rownames(matDistY_sym)

  sharedPopulations <- populationsX[which(populationsX %in% populationsY)]
  sharedPopulations <- sharedPopulations[which(sharedPopulations %in% populationsX)]

  matDistX_pop <- matDistX_sym[sharedPopulations, sharedPopulations]
  matDistY_pop <- matDistY_sym[sharedPopulations, sharedPopulations]


  ### Mantel test ####

  mantelTest <- mantel(matDistX_pop, matDistY_pop)
  mantel_cor <- mantelTest$statistic
  pvalue <- mantelTest$signif

  # Significance code

  if(pvalue >= 0.05){
    signif <- ""
  }
  if(pvalue < 0.05){
    signif <- "."
  }
  if(pvalue < 0.01){
    signif <- "**"
  }
  if(pvalue < 0.001){
    signif <- "***"
  }

  ### Plot ####

  # Scatterplot

  scatter_dist(mat_ld = matDistX_pop,
               mat_gd = matDistY_pop,
               method = "lm",
               thr_gd = NULL,
               thr_ld = NULL,
               se = TRUE,
               smooth_col = "black",
               pts_col = "#999999") +
    xlab(matDistXLab) +
    ylab(matDistYLab) +
    theme_minimal() +
    annotate(geom="text",
             x= quantile(matDistX_sym, probs = 0.5),
             y=max(matDistY_sym),
             label=paste0("r = ",
                          round(mantel_cor, 3), signif),
             color="black")

}
