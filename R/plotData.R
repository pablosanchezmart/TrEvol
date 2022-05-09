#' Plot data on the phylogeny
#'
#' @param phylogeny (phylo) Phylogeny with tip labels contained in dataset$animal
#' @param dataset (data frame) Dataset containing the variable of interest and a column named animal describing terminal taxa of phylogeny.
#' @param variables (character) Names of the variables. They must be contained in dataset.
#' @param colorBy (character) name of the variable to colour the tree by.
#' @param nameBy (character) name of the taxonomic variable to name the tree by.
#' @param panelSpace (numeric)
#' @param FontsizeFactor (numeric)
#'
#' @return
#' @export
#'
#' @examples
plotData <- function(phylogeny, dataset, variables, colorBy = NULL, nameBy = NULL, panelSpace = 4, FontsizeFactor = 0.1){

  if(!requireNamespace("ggtree",  quietly = TRUE)){
    BiocManager::install("ggtree", update = F, ask = F)
  }

  dataset <- dataset[stats::complete.cases(dataset[, variables]), ]
  phylogeny <- ape::drop.tip(phylogeny, phylogeny$tip.label[!phylogeny$tip.label %in% dataset$animal])
  dataset <- as.data.frame(dataset)
  row.names(dataset) <- dataset$genus
  dataset <- dataset[phylogeny$tip.label, ]
  dataset <- dataset[, c("genus", names(dataset))] # we need it as a first column


  p <- ggtree::ggtree(phylogeny)

  if(!is.null(nameBy)){
    taxonList <- unique(dataset[, nameBy])
    taxonList <- taxonList[which(!is.na(taxonList))]

    # name orders
    for(o in taxonList){
      spp <- rownames(dataset[which(dataset[, nameBy] == o), ])
      if(length(spp) == 1){                                               # Some orders have only one species, so MRCA cannot be found
        nd <- which(phylogeny$tip.label == spp)
      } else{
        nd <- phytools::findMRCA(phylogeny, spp, type = "node")
      }
      p <- p + ggtree::geom_cladelabel(node= nd, label= o, fontsize = length(rownames(dataset[dataset$order == o, ])) * FontsizeFactor + 1) +
        ggplot2::coord_cartesian(clip = 'off')  + ggtree::theme_tree2(plot.margin= ggtree::margin(6, 6, 6, 6))
    }
  }

  if(!is.null(colorBy)){
    cat.df <- data.frame(as.factor(dataset[, colorBy]))

    for(var in variables){
      d <- data.frame(id = row.names(dataset), varX = dataset[, var], order = cat.df[, 1])
      p <- ggtree::facet_plot(p, panel = var, data = d, geom = ggstance::geom_barh,
                      mapping =  aes(x = varX, fill = order),
                      stat = "identity")
    }
  } else {
    for(var in variables){
      d <- data.frame(id = row.names(dataset), varX = dataset[, var], order = cat.df[, 1])
      p <- ggtree::facet_plot(p, panel = var, data = d, geom = ggstance::geom_barh,
                      mapping =  aes(x = varX, fill = order),
                      stat = "identity")
    }
  }

  p <- p  +  ggplot2::scale_fill_viridis_d() +
    ggplot2::theme(panel.spacing = unit(panelSpace, "lines")) +  ggplot2::theme(legend.position = "none")
  print(p)
}