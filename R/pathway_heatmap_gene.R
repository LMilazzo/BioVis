# GENE HEATMAP!
#
# This is a function, it filters data a vector of searched genes and returns plots.
#'
#' This function filters a dataset based on provided genes and returns plots
#'
#' @description This function allows you to filter a dataset of pathway enrichment data by a character vector of genes, a ggplot and plotly are both returned respectively in a list.
#'
#' @param data Clustered results from pathfindR experiment using pathfindR cluster_enriched_terms(data, plot_dend = TRUE, plot_clusters_graph = TRUE)
#' @param genes A vector of gene names to filter the data. Default is NULL.
#' @return A list (ggplot, plotly) displaying the data
#' @import dplyr
#' @import ggplot2
#' @import pathfindR
#' @import purrr
#' @import plotly
#' @import htmlwidgets
#' @examples
#' # Example usage:
#' # geneheatmap(data, genes = c("Gene1", "Gene2"))
#' @export
geneheatmap <- function(

    data = NULL,
    genes = NULL

    ){

  #______________INPUT VALIDATION_____________
  #----

  if( is.null(data) ){
    stop("no data input")
  }

  #Ensure required format is met
  expected_names <- c("ID", "Term_Description", "Fold_Enrichment", "occurrence", "support",
                      "lowest_p", "highest_p", "non_Signif_Snw_Genes", "Up_regulated",
                      "Down_regulated", "all_pathway_genes", "num_genes_in_path", "Cluster", "Status")

  if( length(intersect(expected_names, colnames(data))) < length(expected_names)) {
    stop("The data is not the correct format")
  }


  #Check is a list of genes was passed
  if( is.null(genes) ){
    stop("try using a list of genes to plot")
  }

  #----

  #______________DATA MANIPULATION____________
  #----

  #Change the all pathway gene lists generated by pathfindR to a list of individual gene names
  data <- data %>%
    mutate(gene_list = lapply(all_pathway_genes,
                              function(x) trimws(unlist(strsplit(x, ",")))))

  #Initialize a data frame that will hold the plotable information
  tofilter <- data %>% filter(ID == "INITIALIZE")

  #Select any rows of data where any genes in the seached gene list appear
  tofilter <- data %>% filter(map_lgl(gene_list, ~ any(.x %in% genes)))

  #To plot the data there must be at least 2 pathways, if 2 pathways are not present additional data will be
  #added until there is 10 rows
  if(nrow(tofilter) < 2){
    print("your genes did not appear in enough pathways, plot will be populated with top ten pathways")
    data <- data %>% arrange(desc(round(Fold_Enrichment, digits=2)), lowest_p)
    filtered <- rbind(tofilter, head(data, 10)) %>% distinct()
  }else{
    filtered <- tofilter
  }

  #Remove the no longer needed column
  filtered <- filtered %>% select(-gene_list)

  #Reduce to a max of 50 pathways for readability
  if(nrow(filtered) > 50 ){
    print("first 50 pathways will be displayed")
    num_Terms <- 50
  }else{
    num_Terms <- nrow(filtered)
  }

  #Use the pathfindR function to retrieve plot data to make our own plot
  dummy <- term_gene_heatmap(filtered, num_terms = num_Terms, use_description = TRUE)
  data <- dummy$plot_env$term_genes_df
  rm(dummy)

  #Filtered to searched genes
  data <- data %>% filter(Symbol %in% genes)
  #Remove bad rows
  data <- data %>% drop_na()

  if(nrow(data) < 1){
    stop("not enough data to plot")
  }

  w <- nrow(data %>% distinct(Symbol)) * 85

  #----

  #______________CREATE THE PLOT______________
  #----
  plot <- ggplot(data, aes(x = Symbol, y = Enriched_Term, fill = value)) +
    geom_tile(color = 'black') +
    scale_fill_manual(values = c("up" = 'skyblue', "down" = 'indianred3')) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank())

  plotly_plot <- ggplotly(plot, height = 800, width = w) %>%
    layout(
      yaxis = list(
        automargin = TRUE,
        showgrid = TRUE,        # Ensure horizontal grid lines are shown
        gridcolor = "lightgray", # Set grid line color
        gridwidth = 1)
    )
  #----

  return(list(plot, plotly_plot))

}
