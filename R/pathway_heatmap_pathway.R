# PATHWAY HEATMAP!
#
# This is a function, it filters data a vector of searched pathways and returns plots.
#'
#' This function filters a dataset based on provided pathways and returns plots
#'
#' @description This function allows you to filter a dataset of pathway enrichment data by a character vector of pathways or a desired number of pathways to visuallize, a ggplot and plotly are both returned respectively in a list.
#'
#' @param data Clustered results from pathfindR experiment using pathfindR cluster_enriched_terms(data, plot_dend = TRUE, plot_clusters_graph = TRUE)
#' @param pathways Either a numeric that determines the number of genes to visualize or a char vector of gene names to visualize. Default is NULL.
#' @return A list (ggplot, plotly) displaying the data
#' @import dplyr
#' @import ggplot2
#' @import pathfindR
#' @import purrr
#' @import plotly
#' @import htmlwidgets
#' @examples
#' # Example usage:
#' # pathwayheatmap(data, pathways = c("pathway1", "pathway2"))
#' # pathwayheatmap(data, pathways = 20)
#' @export
pathwayheatmap <- function(

    data = NULL,
    pathways = NULL

    ){

  #______________INPUT VALIDATION_____________
  #----

  if( is.null(data) ){
    stop("no input data")
  }

  #Ensure required format is met
  expected_names <- c("ID", "Term_Description", "Fold_Enrichment", "occurrence", "support",
                      "lowest_p", "highest_p", "non_Signif_Snw_Genes", "Up_regulated",
                      "Down_regulated", "all_pathway_genes", "num_genes_in_path", "Cluster", "Status")

  if( length(intersect(expected_names, colnames(data))) < length(expected_names)) {
    stop("The data is not the correct format")
  }

  #----

  #______________DATA MANIPULATION____________
  #----

  #Logic for if nothing is passed or a character vector of pathways is passed
  if( is.null(pathways) || is.character(pathways)){

    if( is.null(pathways) || length(pathways) < 2){
      print("not enough pathways in input, will use top ten")
    }

    #Select the searched paths
    selected_pathways <- data %>%
      filter(Term_Description %in% pathways)

    #select the top ten
    topten <- head(data, 10)

    #if less than 2 are searched then the top ten will always appear
    if(nrow(selected_pathways) < 2){
      data <- rbind(selected_pathways, topten)
    }
    #if they searched enough just those are used
    else{
      data <- selected_pathways
    }
  }

  #If a number is passed in the pathway parameter then the number of pathways use is = to that number up to 75 max
  else if( is.numeric(pathways) ){
    if(pathways > 75){
      print("rounded down to 75")
      pathways <- 75
    }
    data <- head(data, pathways)

  }

  #Stop if the pathway parameter was no good
  else{
    stop("the data input was an unreadable format")
  }

  #just the number of pathways in the df
  terms_to_plot <- nrow(data)

  #Create a dataset using the pathfindR heatmap function
  dummy <- term_gene_heatmap(data, num_terms = terms_to_plot, use_description = TRUE)
  data <- dummy$plot_env$term_genes_df
  dummy <- NULL

  #Drop NAs to speed up plotly
  data <- data %>% drop_na()

  #depreciated feature used to facet by alphabetical group
  # data <- data %>%
  #   mutate(alpha = factor(
  #     case_when(
  #       toupper(substr(Symbol, 1, 1)) %in% LETTERS[1:5] ~ "A-E",
  #       toupper(substr(Symbol, 1, 1)) %in% LETTERS[6:10] ~ "F-J",
  #       toupper(substr(Symbol, 1, 1)) %in% LETTERS[11:15] ~ "K-O",
  #       toupper(substr(Symbol, 1, 1)) %in% LETTERS[16:20] ~ "P-T",
  #       toupper(substr(Symbol, 1, 1)) %in% LETTERS[21:26] ~ "U-Z"
  #     ),
  #     levels = c("A-E", "F-J", "K-O", "P-T", "U-Z")
  #   ))

  #recator the gene column so that the factors are in alphabetical order
  data$Symbol <- factor(data$Symbol, levels = sort(levels(data$Symbol), decreasing=TRUE))

  #dynamically define height of the plot so that each gene gets a min of 20 pixels added to the verticle height
  h <- nlevels(data$Symbol) * 20
  w <- nlevels(data$Term_Description) * 25

  #______________CREATE THE PLOT______________
  #----

  plot <- ggplot(data, aes(x = Enriched_Term, y = Symbol, fill = value)) +
    geom_tile(color = 'black', linewidth = 0.25) +
    scale_fill_manual(values = c("up" = 'skyblue', "down" = 'indianred3')) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 9, angle = 90, margin = margin(0,5,0,5, "pt")),
          axis.title.x = element_blank())


  plotly_plot <- ggplotly(plot, width = w, height = h, dynamicTicks = FALSE) %>%
    layout(
      margin = list(l = 50, r = 50, b = 100, t = 50), # Left, Right, Bottom, Top margins
      yaxis = list(
        automargin = TRUE,
        showgrid = TRUE,        # Ensure horizontal grid lines are shown
        gridcolor = "lightgray", # Set grid line color
        gridwidth = 1    ), # Ensure y-axis labels fit
      xaxis = list(
        showline = FALSE,   # Set x-axis title
        side = "top",              # Manually place title on top
        anchor = "y",              # Anchor x-axis to y-axis
        overlaying = "x",          # Overlay axes (ensures only the top one shows)
        showticklabels = TRUE,
        tickfont = list(size = 12),
        title = ""
      )
    )

  #----

  return(list(plot, plotly_plot))

}
