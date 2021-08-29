#' @title showColorsWithNames
#'
#' @description A quick and dirty way to show colours with name labels in a plot, slightly modified from show_col() function in `scales` library .
#'
#' @param colours A character vector of colours, with or without names
#' @param label_names Whether to show label names
#' @param borders Border color for each tile. Default uses par("fg"). Use border = NA to omit borders.
#' @param cex_label Size of printed labels, as multiplier of default size.
#'
#' @export
#'
showColorsWithNames <- function(colours, label_names = TRUE, cex_label = 1, borders = NULL)
{
  if(is.null(names(colours))){label_names = FALSE}
  n <- length(colours)
  ncol <- ceiling(sqrt(length(colours)))
  nrow <- ceiling(n/ncol)
  colours_plot <- c(colours, rep(NA, nrow * ncol - length(colours)))
  names_label <- c(names(colours), rep("", nrow * ncol - length(colours)))
  colours_plot <- matrix(colours_plot, ncol = ncol, byrow = TRUE)
  names_label <- matrix(names_label, ncol = ncol, byrow = TRUE)
  old <- par(pty = "s", mar = c(0, 0, 0, 0))
  on.exit(par(old))
  size <- max(dim(colours_plot))
  plot(c(0, size), c(0, -size), type = "n", xlab = "", ylab = "",
       axes = FALSE)
  rect(col(colours_plot) - 1, -row(colours_plot) + 1, col(colours_plot), -row(colours_plot),
       col = colours_plot, border = borders)

  hcl <- farver::decode_colour(colours_plot, "rgb", "hcl")
  label_col <- ifelse(hcl[, "l"] > 50, "black", "white")
  if(label_names){
    text(col(colours_plot) - 0.5, -row(colours_plot) + 0.5, paste(colours_plot, names_label,sep = "\n"),
         cex = cex_label, col = label_col)
  } else {
    text(col(colours_plot) - 0.5, -row(colours_plot) + 0.5, colours_plot,
         cex = cex_label, col = label_col)
  }
}
