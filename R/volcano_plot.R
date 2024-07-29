# VOLCANO!
#
# This is a FUNction, it makes a volcano plot!!
# It is a ggplot object.
#
#' Create a Volcano Plot for Gene Expression
#'
#' This function generates a volcano plot to visualize gene expression data.
#'
#' @param data A data frame containing gene expression data.
#' @param down_reg Numeric, threshold for down-regulation (default: -1).
#' @param up_reg Numeric, threshold for up-regulation (default: 1).
#' @param pval Numeric, threshold for p-value significance (default: 0.05).
#' @param population Numeric, proportion of the population to consider (default: 0.3).
#' @param lab_density Numeric, density of labels to be displayed (default: 0.25).
#' @param genes Character vector, specific genes to search and highlight (default: "NOTHING TO SEARCH").
#' @param title Character, title of the plot (default: "").
#' @param subtitle Character, subtitle of the plot (default: "").
#' @param caption Character, caption of the plot (default: "").
#' @return A volcano plot object.
#' @import dplyr
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
erupt <- function(

  #@Param
  data = NULL,
  down_reg = -1,
  up_reg = 1,
  pval = 0.05,
  population = 0.3,
  lab_density = 0.25,
  genes = NULL,
  title = "",
  subtitle = "",
  caption = ""

  ){
#______________INPUT VALIDATION_____________

#DATA
  if( ! is.data.frame( data ) ){
    stop("Data must be a data frame")
  }

  req_col <- c("gene_name", "padj", "log2FoldChange")
  if( length( intersect( colnames(data) , req_col ) ) != 3 ) {
    stop("Missing data in data. data requires 'gene_name', 'padj', 'log2FoldChange'")
  }

  data <- data
  start_rows <- nrow(data)

#NUMERICS
  down_reg <- as.numeric(down_reg)
  up_reg <- as.numeric(up_reg)
  pval <- as.numeric(pval)
  population <- as.numeric(population)
  lab_density <- as.numeric(lab_density)

  if( any( is.na( c(down_reg, up_reg, pval, population, lab_density) ) ) ){
    stop("`down_reg`, `up_reg`, `pval`, `population`, and `lab_density` must be numeric or convertible to numeric.")
  }

#SEARCH GENES
  if( is.null( genes ) ||  length ( genes ) == 0 || all(genes == "") ){
    genes <- character(0)
  }


#______________Handle pvals of 0_____________
  #They are set to 1e-300
  data <- data %>%
    mutate( padj = ifelse( padj <= 1e-300 , 1e-300 , padj ) )

#_____Create a row with expression direction with a character representation____
#down regulated genes are "DOWN"
#up regulated genes are "UP"
#genes found in the search vector are "FOUND"
#genes between expression cutoffs are "NO"

  data <- data %>%

    mutate( ex = case_when(

    gene_name %in% genes ~ "FOUND",

    log2FoldChange > up_reg & padj < pval ~ "UP",

    log2FoldChange < down_reg & padj < pval ~ "DOWN",

    .default = as.character("NO"))

  )


#___Exclude genes found in the search vector from downward filtering___

  found_genes <- NULL

  if(! is.null(genes) ){

    found_genes <- data %>%
      filter(ex == "FOUND") %>%
      mutate(glabel = gene_name)

    data <- data %>% filter(ex != "FOUND")

  }

#___Filtering and sorting the data_____

  data <- data %>%
    arrange(
      ex == "NO", #Genes with no expression are given lowest priority (put at bottom)
      padj,       #Genes with lowest Pvalues have highest priority (put at top)
      desc( abs( log2FoldChange ) ) #fold change is used in case where there is similar pvalue
    )

#____Extract a percentage of the data to plot given by the population argument___
#In case of 0%, 20 genes are plotted anyways

  data <- head( data , ( ( round(nrow(data) * population) ) + 20 ) )

#___Extract a percentage of the genes to name based on the lab density argument___
# 100% will always plot a maximum of 230 labels
# any  genes found withing the top 200 genes by pvalue
# and any genes found within top 30 genes by foldchange
# (with overlap)
# This value can be increased past 100% to label more than 230 genes but plot may
# get cluttered
#
# Genes selected by P value MUTS also have some expression
#
# NOTE: these dataframes only consist of gene names
  p_dens <- round( 10 + (200 - 10) * lab_density )
  l_dens <- round( 5 + (30 - 5) * lab_density )

  p_top <- data %>%
              filter( ex != "NO" ) %>%
              arrange(padj) %>%
              pull(gene_name) %>%
              head(n = p_dens)

  l_top <- data %>%
              arrange( desc(abs(log2FoldChange)) ) %>%
              pull(gene_name) %>%
              head(n = l_dens)

#_____Create a label column in the dataframe___
#Genes found in either p_top or l_top are labeled by their gene_name
#Genes not found in either are NA

  data$glabel <- ifelse(

    data$gene_name %in% c(p_top, l_top),

    data$gene_name, #value for if

    NA #value for else
  )

#The Volcano Plot!!!!!!!!!___________________________________________________

  volcano <- ggplot( data,
                     aes(x = log2FoldChange,
                         y = -log10(padj),
                         color = ex,
                         label = glabel
                        )
                     ) +

    #Plot dot and text labels
    geom_point(size=2) +
    geom_text_repel(max.overlaps = Inf, na.rm=TRUE)+

    #Plot scale
    scale_x_continuous(breaks = seq(-35, 35, 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian() +

    #Color scale
    scale_color_manual(values = c("DOWN" = "#2171b5", "NO" = "grey", "UP" = "#bb0c00"),
                       labels = c("DOWN" = "Downregulated",
                                  "NO" = "Not significant",
                                  "UP" = "Upregulated")) +

    #Verticle lines
    geom_vline(xintercept = c(down_reg, up_reg), color='black', linetype = 'dashed') +

    #Horizontal Pvalue line
    geom_hline(yintercept= -log10(pval), color='black', linetype='dashed') +

    #Axis labels
    labs(x = expression("logFC"),
         y = expression("-log"[10]*"p-value")) +

    #Color Legend label
    labs(color = 'Genes') +

    #Title and subtile
    labs(title = title,
         subtitle = subtitle) +

    #The number of genes displayed of total genes in dataset, And the caption
    labs(caption = paste0("Showing ", nrow(data) + nrow(found_genes), "/", start_rows, " genes",
                          "   --   P-val capped at 1e-300",
                          "\n\n", caption)) +

    theme_minimal() +

    theme(
          plot.margin = margin(10, 10, 10, 10, "pt"),
          axis.title.x = element_text(color='black', size = 20, margin = margin(15, 15, 15, 15, "pt")),
          axis.title.y = element_text(color='black', size=20, margin = margin(15, 15, 15, 15, "pt")),
          axis.text.y = element_text(color='black', size=15),
          axis.text.x = element_text(color='black', size=15),
          axis.line.x = element_line(color='grey'),
          axis.line.y = element_line(color='grey'),

          legend.text = element_text(color='black', size=20),
          legend.title = element_text(color='black', size=20),
          legend.margin = margin(15, 15, 15, 15, "pt"),

          plot.title = element_text(color='black', size=30, margin=margin(20,20,5,10,"pt")),
          plot.subtitle = element_text(color='black', size=20, margin=margin(5,5,15,10,"pt")),
          plot.caption = element_text(color='black', size=15, margin=margin(10, 10, 10, 10, "pt"))
    ) 

  

    #Finally plot genes that were searched on top as a distinct color
  if( ! is.null(found_genes) ){
    volcano <- volcano + geom_point(data = found_genes,
               aes(x = log2FoldChange,
                   y = -log10(padj)
                   ),
               color='lawngreen',
               size=3,
               show.legend = FALSE) +

    geom_text_repel(data = found_genes,
                    aes(x = log2FoldChange,
                        y = -log10(padj),
                        label = gene_name),
                    color='mediumseagreen',
                    show.legend = FALSE)
    }

#Return the ggplot from the function
  return(volcano)

}
