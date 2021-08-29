#' @title Silhouette Plot
#'
#' @description Unsupervised way to choose the optimal clustering resolutions or number of clusters.
#' Computes distance matrix based on correlation distance and calculate silhouette scores for a given clustering result.
#' Outputs a silhouette class object and produces a Silhouette plot.
#'
#' @param object Seurat object.
#' @param cluster_ident Identity class of clustering.
#' @param data_from By default, data_from = "pca", extract pca embedding as input.
#' In case of data_from = "expression", "data" slot of the expression matrix is used.
#'
#'
#' @export
#'
#'
Silhouette <- function(object, cluster_ident, data_from = "pca")
{
  if(data_from == "pca") {
    mat <- t(as.matrix(Seurat::Embeddings(object, reduction = data_from)))
  } else if (data_from == "expression") {
    mat <- as.matrix(Seurat::GetAssayData(object, slot = "data", assay = "RNA"))  # "data" slot of Seurat object
  }
  clust <- object[[cluster_ident]]
  d <- as.dist(1 - cor(mat, method = "pearson"))
  sil=cluster::silhouette(as.numeric(clust[,]), dist = d)
  plot(sil, col = 1:length(levels(clust[,])), border = NA, main = cluster_ident)
  return(sil)
}

#' @title Bootstrapped silhouette score
#'
#' @description Unsupervised way to choose the optimal clustering resolutions or number of clusters.
#' Perform bootstrapping (repeated random subsampling with replacement) to achieve better estimation of Silhouette score.
#'
#' @param object Seurat object.
#' @param cluster_ident Identity class of clustering.
#' @param data_from By default, data_from = "pca", extract pca embedding as input data matrix. In case of data_from = "expression", extract data from original expression matrix.
#' @param balanced Whether or not to perform balanced subsampling, which selects proportional size of cells from each cluster. Default is TRUE.
#' @param iteration Number of iterations for bootstrapping. Default = 100.
#' @param subsample.size Sample size for bootstrapping. If unspecified, subsample one fifth (1/5) of the total population.
#' @param replace Sampling with or without replacement. Default is TRUE.
#' @param verbose Display progress information.
#'
#'
#' @export
#'
#'
SilhouetteBoot <- function(object, cluster_ident, data_from = "pca", balanced = TRUE, iteration = 100, subsample.size = NULL,
                           replace = TRUE, verbose = TRUE)
{
  if(data_from == "pca") {
    mat <- t(as.matrix(Seurat::Embeddings(object, reduction = data_from)))
  } else if (data_from == "expression") {
    mat <- as.matrix(Seurat::GetAssayData(object, slot = "data", assay = "RNA"))  # On expression
  }
  clust <- object[[cluster_ident]]
  out=as.data.frame(matrix(rep(NA, iteration*ncol(clust)), nrow=iteration))
  subsample.size <- switch(is.null(subsample.size)+1, subsample.size, ceiling(ncol(mat)/5))  ## Clever way to parse NULL variables

  if(balanced){
    for(x in 1:iteration){
      if(verbose && x%%10==0){print(paste0("Iteration:",x))}
      for(j in 1:ncol(clust)){
        i=c()
        a=clust[,j]
        for(lab in unique(a)){
          i=c(i,sample((1:ncol(mat))[which(a==lab)],length(which(a==lab))/length(a)*min(subsample.size,ncol(mat))))
        }
        d=as.dist(1 - cor(mat[,i], method = "pearson"))
        if(length(table(clust[i,j]))==1){out[x,j]=0}
        else{
          sil=cluster::silhouette(as.numeric(clust[i,j]),d)
          out[x,j]=mean(sil[, "sil_width"])}}
    }
    means=apply(out,2,median)
    sd=apply(out,2,mad)
    return(list(means,sd))
  } else {
    for(x in 1:iteration)
    {
      i <- unique(sample(1:ncol(mat),min(subsample.size,ncol(mat)), replace = replace))
      d <- as.dist(1 - cor(mat[,i], method = "pearson"))

      for(j in 1:ncol(clust)){
        if(length(table(clust[i,j]))==1){out[x,j]=0}
        else{
          sil=cluster::silhouette(as.numeric(clust[i,j]),d)
          out[x,j]=mean(sil[, "sil_width"])}
      }
    }
    means=apply(out,2,mean)
    sd=apply(out,2,sd)
    return(list(means,sd))
  }
}

