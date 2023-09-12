#' Return complete data for phylogeny and tabular data
#'
#'Return observations present in a data set and also in a phylogeny. Taxon name (e.g., species names) need to be present in the data frame provided.
#'The data frame can contain multiple observations per taxon.
#'
#' @param phylogeny (*phylo*). Phylogeny with tip labels naming terminal taxa.
#' @param dataset (*dataframe*). Data frame with a column representing the terminal taxon names equivalent to tip labels in phylogeny (e.g., species names if working at the species level)
#' @param traits (*character*). Name of the variables of interest. Complete data for these columns will be returned. It is important
#' that the traits specified match column names in the dataset.
#'
#' @return List containing a phylogeny and a data frame including the variables of interest with observations for these variables and
#' also present in the phylogeny.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get phylogenetic and traits complete data for an example obtained simulating data
#'
#' # Simulate example data
#' simulated_traits.data <- simulateDataSet()
#'
#' # Get complete data and phylogeny
#' complete.data <- completePhyloData(
#' phylogeny = simulated_traits.data$phylogeny,
#' dataset = simulated_traits.data$data,
#' traits = colnames(simulated_traits.data$data),
#' taxon_name = "animal")
#'
#' # Data frame with complete data for the specified traits
#' complete.data$dta
#'
#' # Phylogeny for the complete data
#'
#' complete.data$phylo
#' }
completePhyloData <- function(phylogeny, dataset, traits, taxon_name = "animal") {
  phyData <- list()

  # Make sure dataset is in data.frame format
  dataset <- as.data.frame(dataset)

  # Select complete cases for the variables specified in the argument traits
  completeVec <- stats::complete.cases(dataset[, traits])
  completeData <- dataset[completeVec, ]

  # Get a phylogeny for the complete data
  phylo <- ape::keep.tip(phylogeny, phylogeny$tip.label[phylogeny$tip.label %in% completeData[, taxon_name]])
  phyData$dta <- as.data.frame(completeData[completeData[, taxon_name] %in% phylo$tip.label, ])

  # Ordered taxa as in the phylogeny
  taxon_ordered <- data.frame(n = 1:length(phylo$tip.label))
  taxon_ordered[, taxon_name] <- phylo$tip.label
  taxon_ordered$n <- NULL

  # Make sure that tip nodes are deleted as they can give some problems
  if(!is.null(phylo$tip.nodes)){
    phylo$tip.nodes <- NULL
  }

  # rename family as it gives problems with some other functions of the packages coming from MCMCglmm
  if("family" %in% colnames(phyData$dta)){
    phyData$dta <- phyData$dta %>% dplyr::rename(Family = family)
    print("renaming family column to Family to avoid problems with MCMCglmm")
  }

  # Join data with ordered taxa (ordered as in the phylogeny)
  phyData$dta <- as.data.frame(dplyr::left_join(taxon_ordered, phyData$dta, by= taxon_name))
  phyData$phylo <- phylo

  return(phyData)
}
