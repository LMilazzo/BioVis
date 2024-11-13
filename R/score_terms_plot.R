# SCORE TERMS!
#
# This is a function, it filters data and plots the score of pathway terms facet by case vs control.
#'
#' This function filters a dataset based on provided pathways
#'
#' @description This function allows you to filter a dataset of pathway enrichment data by a character vector of pathways, a ggplot and recommended dimensions are returned in a list.
#'
#' @param data Clustered results from pathfindR experiment using pathfindR cluster_enriched_terms(data, plot_dend = TRUE, plot_clusters_graph = TRUE)
#' @param pathways A vector of pathway names to filter the data. Default is NULL.
#' @param abundance a dataframe of normalized counts must contain 1 column of Gene_symbol
#' @param repOnly True if to only plot the representative status pathways
#' @param cases A vector of sample names that belong to the experiment case
#' @return A list (ggplot, height, width) displaying the data
#' @import dplyr
#' @import ggplot2
#' @import pathfindR
#' @examples
#' # Example usage:
#' # score_pathways(data, abundance, cases, repOnly = TRUE, pathways = c('pathway1'))
#' @export
score_pathway_terms <- function(

  data = NULL,
  abundance = NULL,
  cases = NULL,
  repOnly = TRUE,
  pathways = NULL

){

  if(is.null(data)){
    stop('No data was passed')
  }

  if(is.null(abundance)){
    stop('No abundance data')
  }

  #Check data types
  original_na_counts <- sapply(colnames(data),
                               function(col) sum(is.na(data[[col]])))

  data$ID <- as.character(data$ID)
  data$Term_Description <- as.character(data$Term_Description)
  data$non_Signif_Snw_Genes <- as.character(data$non_Signif_Snw_Genes)
  data$Up_regulated <- as.character(data$Up_regulated)
  data$Down_regulated <- as.character(data$Down_regulated)
  data$all_pathway_genes <- as.character(data$all_pathway_genes)
  data$Status <- as.character(data$Status)
  data$Fold_Enrichment <- as.numeric(data$Fold_Enrichment)
  data$occurrence <- as.numeric(data$occurrence)
  data$support <- as.numeric(data$support)
  data$lowest_p <- as.numeric(data$lowest_p)
  data$highest_p <- as.numeric(data$highest_p)
  data$num_genes_in_path <- as.numeric(data$num_genes_in_path)
  data$Cluster <- as.numeric(data$Cluster)

  for(i in colnames(data)){
    new_na <- sum(is.na(data[[i]]))
    if(!new_na == original_na_counts[[i]]){
      stop(paste0('Coercion of column ', i, ' to proper class introduced na values. ', i, ' requires numeric.'))
    }
  }

  #Check abundance
  if(! 'Gene_symbol' %in% colnames(abundance)){
    stop('Gene_symbol is a req column in abundance data')
  }
  print(head(abundance))
  rownames(abundance) <- abundance$Gene_symbol
  abundance <- abundance %>% select(-Gene_symbol)
  abundance <- as.matrix(abundance)

  print(head(abundance))
  if(!is.numeric(abundance)){
    stop('non numeric value found in abundance data matrix')
  }

  #_________________________________________________________________________#

  #filter if
  if(repOnly == TRUE){
    data <- data %>% filter(Status == 'Representative')
  }

  #Make scores
  scores <- score_terms(
    enrichment_table = data,
    exp_mat = abundance,
    cases=cases,
    use_description = TRUE, # default FALSE
    label_samples = TRUE, # default = TRUE,
    plot_hmap = FALSE)

  #more filter for pathways
  if(!is.null(pathways)){
    scores <- data.frame(scores) %>% mutate(term = rownames(scores))
    scores <- as.matrix(scores[scores$term %in% pathways,  , drop = FALSE] %>% select(-term))
  }

  temp <- data.frame(scores)
  h <- nrow(temp) * 20
  if(h < 600){
    h <- 600
  }
  w <- ncol(temp) * 50
  if(w < 600){
    w <- w + 600
  }

  print(scores)

  #make plot
  plot <- plot_scores(
    score_matrix = scores,
    cases = cases,
    label_samples = TRUE, # default = TRUE
    low = "#f7797d", # default = "green"
    mid = "#fffde4", # default = "black"
    high = "#1f4037" # default = "red"
  ) +
    theme(
      axis.text.x = element_text(size = 15, color = 'black'),
      axis.text.y = element_text(size = 15, color = 'black'),
      panel.background = element_blank()
    ) +
    labs(x = '', y = '')

  return(list(plot, h, w))

}
