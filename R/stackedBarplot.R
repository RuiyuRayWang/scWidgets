#' @title Stacked Barplot
#'
#' @description Plots stacked barplots for single cell datasets
#'
#' @details
#' Draws barplots of single cell data (gene expression, metrics, PC scores, etc.), value of each cell is represented by one bar.
#'
#' @param object A Seurat object.
#' @param features Features to plot (gene expression, metrics, PC scores, anything that can be retreived by FetchData)
#' @param assay Name of assay to use, defaults to "RNA".
#' @param slot Slot to draw data from, defaults to "data"
#' @param group.by Group (facet and color) cells in different ways (for example, orig.ident)
#'
#' @examples
#' data("pbmc_small)
#' StackedBarplot(pbmc_small, features = c("RGS1", "TYROBP", "GZMA", "CCL4", "CD79A", "MS4A1"), group.by = "groups")
#'
#' @return A patchworked ggplot object
#' @export
#'
StackedBarplot <- function(object, features, assay = "RNA", slot = "data", group.by = NULL){
  if(is.null(assay)){assay = Seurat::DefaultAssay(object)}

  df <- SeuratObject::FetchData(object, vars = features, slot = slot)
  df <- tibble::rownames_to_column(df, var = "cells")
  if(!is.null(group.by)){
    groups <- object[[group.by]]
  } else {
    groups <- data.frame(group = rep("", ncol(object)), row.names = colnames(object))
  }
  groups <- tibble::rownames_to_column(groups, var = "cells")

  df <- dplyr::left_join(df, groups, by = "cells")
  df <- tibble::column_to_rownames(df, var = "cells")

  df_long <- tibble::rownames_to_column(df, var = "cells")
  df_long <- tidyr::pivot_longer(df_long, cols = -c("cells",all_of(group.by)), names_to = "features", values_to = "value")

  df_long$value <- as.numeric(df_long$value)
  df_long$features <- factor(df_long$features, levels = features)

  p <- ggplot2::ggplot(df_long, ggplot2::aes_string(x = "cells", y = "value", fill = group.by, color = group.by)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_grid(reformulate(group.by, "features")) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0),
      axis.ticks = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(),
      legend.position = "none"
      )

  return(p)
}
