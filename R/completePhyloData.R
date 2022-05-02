#' Complete phylogeny and traits data
#'
#' @param phylogeny phylogeny with tip labels naming terminal taxa
#' @param dataset dataframe with a column including taxa named "animal" equivalent to tip labels in phylogeny
#' @param traits character vector
#'
#' @return
#' List with phylogeny and traits data including taxa with both information a
#' @export
#'
#' @examples
completePhyloFun <- function(phylogeny, dataset, traits) {
  phyData <- list()
  completeVec <- stats::complete.cases(dataset[, traits])
  completeData <- dataset[completeVec, ]

  phylo <- ape::keep.tip(phylogeny, phylogeny$tip.label[phylogeny$tip.label %in% completeData$animal])
  phyData$dta <- as.data.frame(completeData[completeData$animal %in% phylo$tip.label, ])
  rownames(phyData$dta) <- phyData$dta$animal
  phyData$dta <- phyData$dta[phylo$tip.label, ]
  phyData$phylo <- phylo

  return(phyData)
}
