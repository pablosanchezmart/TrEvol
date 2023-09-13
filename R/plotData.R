#' Plot data on the phylogeny
#'
#'Plot variables as barplots for each tip of the phylogeny with data available. Several variables can be plotted at the same time.
#'
#' @param variables (*character*). Names of the variable or variables to be ploted as they appear in data.
#' @param dataset (*data frame*). Data containing the variable of interest and a column naming the terminal taxa as it appears in the phylogeny (e.g., species names).
#' @param phylogeny (*phylo*). Phylogeny with tip labels contained in data.
#' @param color_by (*character*). Name of the variable to colour the tree by, if desired.
#' @param name_by (*character*). Name of the taxonomic variable to display in the phylogeny panel.
#' @param panel_space (*numeric*). Size of each one of the panels (one for the phylogeny and one per variable included in the "variables" argument).
#' @param fontsize_factor (*numeric*). Font size multiplier. It modifies the size of the font.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate example data
#' simulated_traits.data <- simulateDataSet()
#'
#' # Plot G1_trait1 data on simulated phylogeny
#'
#' plotData(variables = "G1_trait1",
#' dataset = simulated_traits.data$data,
#' phylogeny = simulatedTraits.data$phylogeny,,
#' terminal.taxon = "species"
#' )
#' }
plotData <- function (variables = NULL,
                      dataset = NULL,
                      phylogeny = NULL,
                      terminal.taxon = "species",
                      color_by = NULL,
                      name_by = NULL,
                      panel_space = 4,
                      fontsize_factor = 0.1) {

  # Arguments
  if(is.null(variables)){
    stop("Specify variables argument")
  }

  if(is.null(dataset)){
    stop("Specify dataset argument")
  }


  # Install ggtree if required
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    BiocManager::install("ggtree", update = F, ask = F)
  }

  # Give common name to taxonomy (animal)
  dataset$animal <- dataset[, terminal.taxon]

  # Complete cases
  dataset <- dataset[stats::complete.cases(dataset[, variables]), ]
  phylogeny <- ape::drop.tip(phylogeny, phylogeny$tip.label[!phylogeny$tip.label %in% dataset$animal])


  dataset <- as.data.frame(dataset)
  row.names(dataset) <- dataset$animal
  dataset <- dataset[phylogeny$tip.label, ]

  # plot
  p <- ggtree::ggtree(phylogeny) +  ggplot2::coord_cartesian(clip = "off")

  # Name terminal taxons
  if (!is.null(name_by)) {
    taxonList <- unique(dataset[, name_by])
    taxonList <- taxonList[which(!is.na(taxonList))]
    for (o in taxonList) {
      spp <- rownames(dataset[which(dataset[, name_by] == o), ])
      if (length(spp) == 1) {
        nd <- which(phylogeny$tip.label == spp)
      }
      else {
        nd <- phytools::findMRCA(phylogeny, spp, type = "node")
      }

      p <- p + ggtree::geom_cladelabel(node = nd, label = o, fontsize = length(rownames(dataset[dataset[, name_by] == o, ])) * fontsize_factor + 1) +
        ggtree::theme_tree2(plot.margin = ggtree::margin(6, 6, 6, 6))
    }
  }

  # Colour tree
  if (!is.null(color_by)) {
    cat.df <- data.frame(as.factor(dataset[, color_by]))
    for (var in variables) {
      d <- data.frame(id = row.names(dataset), varX = dataset[, var], color = cat.df[, 1])
      p <- ggtree::facet_plot(p, panel = var, data = d, geom = ggstance::geom_barh, mapping = ggtree::aes(x = varX, fill = color), stat = "identity")
    }
  } else {
    for (var in variables) {
      d <- data.frame(id = row.names(dataset), varX = dataset[, var])
      p <- ggtree::facet_plot(p, panel = var, data = d, geom = ggstance::geom_barh, mapping = ggtree::aes(x = varX), stat = "identity")
    }
  }

  # plot theme
  p <- p + ggplot2::scale_fill_viridis_d() + ggplot2::theme(panel.spacing = ggtree::unit(panel_space, "lines")) + ggplot2::theme(legend.position = "none")

  # plot
  print(p)

  return(p)
}
