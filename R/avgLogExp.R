#' @title Compute average log expression
#'
#' @description Custom function to compute log average expression level
#'
#' @param object
#' @param features
#' @param cells
#' @param pseudocount.use
#' @return Mean in log space
#' @export
#'
avgLogExp <- function(object, features = NULL, cells = NULL, pseudocount.use = 1) {
  if(is.null(cells)){
    cells = Seurat::Cells(object)  # Default is to compute on all cells
  }

  x = Seurat::GetAssayData(object, slot = "data", assay = "RNA")[features,cells]

  if(is.numeric(x)){

    log_mean_exp = log(x = mean(x = expm1(x = x)) + pseudocount.use)

  } else if (is(x, "sparseMatrix")){

  log_mean_exp = log(x = Matrix::rowMeans(x = expm1(x = x)) + pseudocount.use)

  }

  return(log_mean_exp)
}
