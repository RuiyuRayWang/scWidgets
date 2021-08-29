#' @title Geometric Sketching
#'
#' @description Geometric Sketching to subsample data from single cell genomics into smaller representative subset of cells, while preserving
#' biological complexity, highlighting rare cell states, and accelerating complex analysis such as data integration.
#'
#' @details
#' The original package `geosketch` is written in Python. With `reticulate`, we implement it in R and provides a function that directly works with
#' the Seurat object.
#'
#' Publication: Hie, B., Cho, H., DeMeo, B., Bryson, B., and Berger, B. (2019). Geometric Sketching Compactly Summarizes the Single-Cell
#' Transcriptomic Landscape. Cell Systems.
#'
#' Github repo: https://github.com/brianhie/geosketch
#'
#' @param object A Seurat object.
#' @param sketch.size Number of cells to sketch. By default sketch 1000 cells.
#' @param assay Assay used for sketching. "RNA" by default.
#' @param slot By default slot = "data".
#' @param k Target rank of the low-rank SVD. k should satisfy k << min(m,n). By default k = 20.
#'
#' @examples
#' ## If haven't installed already.
#' reticulate::py_install("geosketch")
#'
#' ## mgc has 7527 cells.
#' mgc.sketched <- geosketch(mgc)
#'
#' @return A sketched Seurat object.
#' @export
#'
geosketch <- function(object, sketch.size = 1000, assay = "RNA", slot = "data", k = 20){
  geosketch.py <- reticulate::import('geosketch')
  X <- t(as.matrix(Seurat::GetAssayData(object, assay = assay, slot = slot)))
  # Get top PCs from randomized SVD
  s <- rsvd::rsvd(X, k = k)
  X.pcs <- s$u %*% diag(s$d)

  sketch.size <- as.integer(sketch.size)
  sketch.indices <- geosketch.py$gs(X.pcs, sketch.size, one_indexed = TRUE)

  object.sketch <- subset(object, cells = Seurat::Cells(object)[unlist(sketch.indices)])
  return(object.sketch)
}
