# GENE COUNTS!
#
# This is a function, it plots gene counts!!
# It is a ggplot object.
#
#' Plot Gene Counts for a Single Gene
#'
#' This function generates a plot to visualize normalized counts for a single gene across samples.
#'
#' @description This function takes a single row of gene data and associated metadata to create a plot of normalized counts for each sample.
#'
#' @param gene A data frame containing a single row with gene information including gene name, gene id, fold change, and adjusted p-value.
#' @param metadata A data frame containing metadata for the samples.
#' @param conditions A vector of column names in the metadata to be used for the plot aesthetics.
#' @return A ggplot object representing the gene counts across samples with additional conditions used for coloring, shaping, and sizing the points.
#' @import ggplot2
#' @import dplyr
#' @import beeswarm
#' @examples
#' # Example usage:
#' # gene <- data.frame(gene_name="Gene1", gene_id="G1", padj=0.01, log2FoldChange=2, sample1=10, sample2=15, sample3=5)
#' # metadata <- data.frame(sample=c("sample1", "sample2", "sample3"), condition=c("A", "B", "A"))
#' # plotGeneCounts(gene, metadata, c("condition"))
#' @export
gplop <- function(

  gene = NULL,
  metadata = NULL,
  conditions = NULL,
  textSizeAdjustment = 0

  ){

  #______________INPUT VALIDATION_____________
  #----
  # Check if metadata is provided
  if( is.null( metadata ) ){
    stop("Metadata must be provided")
  }

  # Check if gene is a data frame
  if( !is.data.frame(gene) || is.null(gene) ){
    stop("Gene must be a provided data frame")
  }

  #check if the required data is present
  req_col <- c("gene_name", "padj", "log2FoldChange")
  if( length(intersect(colnames(gene) , req_col )) < length(req_col) ){
    stop("A required data column ('gene_name', 'padj', 'log2FoldChange') is missing")
  }

  if( length(colnames(gene %>% select(starts_with('.')))) < 2){
    stop("Denote sample columns starting with '.', or only one sample")
  }

  #----

  #_________Extract gene information_________
  #----

  name <- ( gene %>% select(gene_name) )[1]

  p <- format( (gene %>% select(padj))[1], scientific = TRUE, digits = 8 )

  fold <- round( (gene %>% select(log2FoldChange))[1] , digits=8 )

  counts <- gene %>% select(starts_with('.'))

  mean <- round( mean( as.numeric(counts) ) )

  #----

  #_________Add counts to metadata_________
  #----

  toplot <- metadata %>% mutate( counts = data.frame(t(counts))[,1] )

  #----

  #_________Identify intersecting conditions_________
  #----

  #if 0 set to the first available cond

  cond <- intersect(conditions, colnames(metadata))
  if( is.null(cond) || length(cond) <= 0 ){
    cond <- metadata[,1]
  }
  #----

  #_________Set up plot aesthetics_________
  #----

  #only the first 4 things are used

  # Set up plot aesthetics
  aesthetics <- aes(x = as.factor(toplot[[cond[1]]]), y = toplot$counts)

  a <- length(cond)
  if (a >= 2) {
    aesthetics <- modifyList(aesthetics, aes(colour = as.factor(toplot[[cond[2]]])))
  }
  if (a >= 3) {
    aesthetics <- modifyList(aesthetics, aes(shape = as.factor(toplot[[cond[3]]])))
  }
  if (a >= 4) {
    aesthetics <- modifyList(aesthetics, aes(size = as.factor(toplot[[cond[4]]])))
  }

  #----

  #_________Create the plot_________
  #----

  plot <- ggplot( toplot, aesthetics ) +
    xlab( cond[1] ) +
    ylab( "Count" ) +

    theme_minimal() +

    labs(title = name,
         subtitle = paste("Fold Change: ", fold),
         caption = paste("Mean Count: ", mean, "   Padj: ", p ),
         color = if (a >= 2) cond[2] else NULL,
         shape = if (a >= 3) cond[3] else NULL,
         size = if (a >= 4) cond[4] else NULL
    ) +

    theme(plot.margin = margin(3, 2, 2, 2, "pt"),
          axis.title.x = element_text(color='black', size = 15 + textSizeAdjustment, margin = margin(2, 2, 2, 2, "pt")),
          axis.title.y = element_text(color='black', size=15 + textSizeAdjustment, margin = margin(2, 2, 2, 2, "pt")),
          axis.text.y = element_text(color='black', size=15 + textSizeAdjustment),
          axis.text.x = element_text(color='black', size=15 + textSizeAdjustment),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line.x = element_line(color='grey'),
          axis.line.y = element_line(color='grey'),
          legend.text = element_text(color='black', size=20 + textSizeAdjustment),
          legend.title = element_text(color='black', size=20 + textSizeAdjustment),
          legend.margin = margin(2, 2, 2, 2, "pt"),
          plot.title = element_text(color='black', size=20 + textSizeAdjustment, margin=margin(2, 2, 2, 2, "pt")),
          plot.subtitle = element_text(color='black', size=15 + textSizeAdjustment, margin=margin(2, 2, 5, 2, "pt")),
          plot.caption = element_text(color='black', size=10 + textSizeAdjustment, margin=margin(2, 80, 2, 1, "pt")),
          plot.margin = margin(8,8,8,8,"pt")
    )

  # Change size of points according to the number of conditions
  if(a >= 4){
    plot <- plot + scale_size_discrete(range = c(3, 10)) + geom_beeswarm(cex = 3)
  }else{
    plot <- plot + geom_beeswarm(size=3, cex = 3)
  }

  #----

  return(plot)

}
