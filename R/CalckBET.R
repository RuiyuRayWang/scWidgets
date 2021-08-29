#' @title k-Nearest Neighbour Batch Effect Test
#'
#' @description Wrapper function of kBET for Seurat objects.
#'
#' @details
#' General principles behind kBET:
#'
#' When batch effect exists in a dataset, the dataset contains disproportional amounts of samples from each batch within neiborhoods surrounding
#' each point. Using Chi-squared statistics, test whether the proportion of each batch within a neighborhood is disproportional. If it's not
#' proportional (i.e. p < critical value), reject null hypothesis (proportional distribution).
#' kBET aggregates the test results computed from multiple neighborhoods, and reports a "rejection rate" as a metric for batch effect. High
#' rejection rate indicates strong batch effect, whereas low "rejection rate" indicates mild batch effect.
#' For "acceptance rate", it is simply a rescaled value of "rejection rate", computed as `acceptance rate = 1 - rejection rate`.
#'
#' Publication: Büttner, M., Miao, Z., Wolf, F.A., Teichmann, S.A., and Theis, F.J. (2019). A test metric for assessing single-cell RNA-seq batch
#' correction. Nat Methods.
#'
#' Github repo: https://github.com/theislab/kBET
#'
#' @param object A Seurat object.
#' @param sketch.size Number of cells to sketch. By default sketch 1000 cells.
#' @param assay Assay used for sketching. "RNA" by default.
#' @param slot By default slot = "data".
#' @param k0 Neighborhood size. By default k0 = mean batch size. To prevent kBET rejection rate from saturating to 1, lower this value.
#' @param ... Arguments passed to kBET.
#'
#' @examples
#' res <- CalckBET(pbmc_small, "groups")
#'
#' @return A list containing detailed results calculated by kBET.
#' @export
#'
CalckBET <- function(object, ident, k0 = NULL, knn = NULL, assay = "RNA", slot = "scale.data",
                      ...){
  batch <- as.character(object[[ident]][,1])
  data <- t(as.matrix(Seurat::GetAssayData(object, assay = assay, slot = slot)))
  batch.estimate <- kBET::kBET(data, batch, k0 = k0, ...)
  return(batch.estimate)
}
