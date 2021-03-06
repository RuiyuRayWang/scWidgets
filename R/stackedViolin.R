#' @title Parse Seurat violin plot
#'
#' @description Helper function to reproduce scanpy stacked violin plot with Seurat.
#' @details
#' With violin plot generated by `Seurat::VlnPlot`, this function removes
#' its x-axis text and tick, and adjusts the white space between each
#' plot with plot.margin argument.
#' ... pass any arguments to VlnPlot in Seurat
#' Credit to: Tang Ming
#' https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#'
#' @param obj A Seurat object.
#' @param features Genes to be plotted.
#' @param pt.size Size of jitter points. Default is 0 (No point is plotted).
#' @param plot.margin Adjust the white space between each plot.
#' @return A ggplot object
#' @export
#'
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = ggplot2::unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- Seurat::VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    ggplot2::xlab("") + ggplot2::ylab(feature) + ggplot2::ggtitle(label = NULL) +
    ggplot2::theme(legend.position = "none",
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), angle = 0),
          axis.text.y = ggplot2::element_text(size = ggplot2::rel(1)),
          plot.margin = plot.margin )
  return(p)
}

#' @title Extract max value of y axis
#'
#' @description  Helper function to reproduce scanpy stacked violin plot with Seurat.
#' @details
#' Credit to: Tang Ming
#' https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#'
#' @param obj A Seurat object.
#' @param features Genes to be plotted.
#' @param pt.size Size of jitter points. Default is 0 (No point is plotted).
#' @param plot.margin Adjust the white space between each plot.
#' @return ymax
#'
## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

#' @title Parse title length
#'
#' @description Helper function to reproduce scanpy stacked violin plot with Seurat.
#' @details
#' Credit to: Tang Ming
#' https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#'
#' @param wrap_title Title to be plotted.
#' @param break_length Length of break point: where to insert line break?
#' @return String of title with appropriate breaks
#'
parse_title <- function(wrap_title, break_length = 23){
  paste(strwrap(wrap_title, break_length), collapse="\n")
}

#' @title Make stacked violin plot
#'
#' @description Reproduce scanpy stacked violin plot with Seurat.
#' @details
#' Credit to: Tang Ming
#' https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#'
#' @param obj A Seurat object.
#' @param features Genes to be plotted.
#' @param pt.size Size of jitter points. Default is 0 (No point is plotted).
#' @param plot.margin Adjust the white space between each plot.
#' @return A ggplot object
#'
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = ggplot2::unit(c(-0.75, 0, -0.75, 0), "cm"),
                          wrap_title = NULL, break_length = 23,
                          ...) {

  if(is.null(wrap_title)){wrap_title = paste0(length(features)," genes")}

  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    ggplot2::theme(axis.text.x=ggplot2::element_text(), axis.ticks.x = ggplot2::element_line())

  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x +
                             ggplot2::scale_y_continuous(breaks = c(y)) +
                             ggplot2::expand_limits(y = y))

  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)

  p <- p + patchwork::plot_annotation(
    title = parse_title(wrap_title, break_length),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold"))
  )
  return(p)
}
