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
completePhyloData <- function(phylogeny, dataset, traits, taxonName = "animal") {
  phyData <- list()

  dataset <- as.data.frame(dataset)

  completeVec <- stats::complete.cases(dataset[, traits])
  completeData <- dataset[completeVec, ]

  phylo <- ape::keep.tip(phylogeny, phylogeny$tip.label[phylogeny$tip.label %in% completeData[, taxonName]])
  phyData$dta <- as.data.frame(completeData[completeData[, taxonName] %in% phylo$tip.label, ])

  taxon_ordered <- data.frame(n = 1:length(phylo$tip.label))
  taxon_ordered[, taxonName] <- phylo$tip.label
  taxon_ordered$n <- NULL

  if(!is.null(phylo$tip.nodes)){
    phylo$tip.nodes <- NULL
  }

  if("family" %in% colnames(phyData$dta)){
    phyData$dta <- phyData$dta %>% dplyr::rename(Family = family)
    print("renaming family column to Family to avoid problems with MCMCglmm")
  }
  phyData$dta <- as.data.frame(dplyr::left_join(taxon_ordered, phyData$dta, by= taxonName))
  phyData$phylo <- phylo

  return(phyData)
}
