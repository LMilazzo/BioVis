# HeatMap Single PATHWAY!
#
# This is a function, it makes a heatmap of a single pathway vs. samples!!
# It is a ggplot object.
#
#' Create a heatmap for Gene Expression
#'
#'
#' @param DEG_results A data frame containing gene expression data.
#' @param Pathfinder_results A data frame containing pathway expression data
#' @param pathway a sting of the searched pathway to match
#' @param genes_listed the number of pathway genes to include ranked by padj
#' @return A ggplot heatmap
#' @import dplyr
#' @import ggplot2
#' @import pheatmap
#' @import ggplotify
#' @export
single_pathway_heatmap <- function(

  pathway = NULL,
  DEG_results = NULL,
  Pathfinder_results = NULL,
  genes_listed = 100

){

  #Write Assertions and data checks

  if(is.null(pathway)){
    stop('Give a pathway to map')
  }

  if(is.null(DEG_results) || is.null(Pathfinder_results)){
    stop('Data needed')
  }

  if(! all(c('padj','gene_name') %in% colnames(DEG_results)) ){
    stop('DEG results needs a padj and gene_name column')
  }

  #Check Pathfinder p columns
  p <- Pathfinder_results
  expected_names <- c("ID", "Term_Description", "Fold_Enrichment", "occurrence", "support",
                      "lowest_p", "highest_p", "non_Signif_Snw_Genes", "Up_regulated",
                      "Down_regulated", "all_pathway_genes", "num_genes_in_path", "Cluster", "Status")

  if( length(intersect(expected_names, colnames(p))) < length(expected_names)) {
    stop("The data is not the correct format")
  }

  #Check p types
  original_na_counts <- sapply(colnames(p),
                               function(col) sum(is.na(p[[col]])))

  p$ID <- as.character(p$ID)
  p$Term_Description <- as.character(p$Term_Description)
  p$non_Signif_Snw_Genes <- as.character(p$non_Signif_Snw_Genes)
  p$Up_regulated <- as.character(p$Up_regulated)
  p$Down_regulated <- as.character(p$Down_regulated)
  p$all_pathway_genes <- as.character(p$all_pathway_genes)
  p$Status <- as.character(p$Status)
  p$Fold_Enrichment <- as.numeric(p$Fold_Enrichment)
  p$occurrence <- as.numeric(p$occurrence)
  p$support <- as.numeric(p$support)
  p$lowest_p <- as.numeric(p$lowest_p)
  p$highest_p <- as.numeric(p$highest_p)
  p$num_genes_in_path <- as.numeric(p$num_genes_in_path)
  p$Cluster <- as.numeric(p$Cluster)

  for(i in colnames(p)){
    new_na <- sum(is.na(p[[i]]))
    if(!new_na == original_na_counts[[i]]){
      stop(paste0('Coercion of column ', i, ' to proper class introduced na values. ', i, ' requires numeric.'))
    }
  }

  #Main stuff

  p$Term_Description <- tolower(p$Term_Description)

  single_p <- p %>%
    filter(Term_Description == tolower(pathway)) %>%
    select(Term_Description, Fold_Enrichment, all_pathway_genes,
           Up_regulated, Down_regulated, Status, non_Signif_Snw_Genes)

  if(nrow(single_p) < 1){
    return()
  }

  if(nrow(single_p) > 1){
    stop('Mulitple matches found')
  }

  up_regs <- strsplit(single_p$Up_regulated, ", ")
  down_regs <- strsplit(single_p$Down_regulated, ", ")
  non_sigs <- strsplit(single_p$non_Signif_Snw_Genes, ", ")

  #make gene list
  gene_list <- data.frame(gene_name = c(unlist(up_regs), unlist(down_regs), unlist(non_sigs)))


  if(genes_listed < 10){
    genes_listed <- 10
  }

  #Sample columns
  data <- DEG_results %>%
    select(padj, gene_name, starts_with('.')) %>%
    filter(tolower(gene_name) %in% tolower(gene_list$gene_name)) %>%
    arrange(padj) %>%
    head(n = genes_listed) %>%
    select(-padj)

  rownames(data) <- data$gene_name
  data <- data %>%select(-gene_name) %>% as.matrix()

  if(ncol(data) < 1){
    stop('No sample columns found')
  }

  data_log2 <- log2(data + 1)

  return(data_log2)

}
