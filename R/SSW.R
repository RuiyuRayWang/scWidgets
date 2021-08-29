#' @title Seurat Standard Workflow (SSW)
#'
#' @description Run Seurat Standard Workflow.
#'
#' - Normalization \cr
#' - Feature Selection (Highly variable genes, HVGs) \cr
#' - Data scaling \cr
#' - PCA \cr
#' - Construct SNN graph \cr
#' - Unsupervised Clustering \cr
#' - Non-linear dimension reduction (tSNE, UMAP) \cr
#'
#' By default, we compute 50 PCs and use the 1~50 PCs as input dimension. To use an optimized dimensionality parameter, run JackStraw or ElbowPlot
#' to determine the dimensionality of the dataset and supply to `dims`.
#'
#' @param object A Seurat object.
#' @param assay "RNA" or "integrated".
#' @param nfeatures Number of HVGs selected.
#' @param PC_features If supplied, run PCA on these features.
#' @param npcs Number of PCs to compute.
#' @param dims Chosen dimensionality of the data. Used in SNN-graph construction and Non-linear dimension reduction.
#' @param k.param Defines k for the k-nearest neighbor algorithm. Supplied to FindNeighbors.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Supplied to Findclusters.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Supplied to Findclusters.
#' @param perplexity Perplexity used in RunTSNE.
#' @param verbose Print output.
#' @examples
#' ## pbmc_small is too small, must set perplexity to lower value
#' pbmc_small <- SSW(pbmc_small, perplexity = 10)
#'
#' @export
#'
SSW <- function(object, assay=NULL, nfeatures = 2000, PC_features=NULL, npcs=50, dims=1:50, k.param=20, algorithm=1, resolution=0.3,
                perplexity = 30, verbose = TRUE) {

  if(is.null(assay)){
    if("integrated" %in% names(object)){
      assay <- "integrated"  # Default is to work on integrated assay
    } else {
      assay <- "RNA"
    }
  }

  if(assay=='integrated'){
    Seurat::DefaultAssay(object) <- 'integrated'
  }else if(assay=='RNA'){
    Seurat::DefaultAssay(object) <- 'RNA'
    object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = verbose)
    object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures, verbose = verbose)
  }
  # Seurat Standard workflow
  object <- Seurat::ScaleData(object, verbose = verbose)
  object <- Seurat::RunPCA(object, features = PC_features, npcs = npcs, verbose = verbose)

  if(assay=='integrated'){
    object <- Seurat::RunUMAP(object, reduction = "pca", dims = dims, umap.method = 'umap-learn', metric = 'correlation', reduction.key = "UMAPint_", reduction.name = "umap.int", verbose = verbose)
    object <- Seurat::RunTSNE(object, reduction = "pca", dims = dims, reduction.key = 'TSNEint_', reduction.name = 'tsne.int', perplexity = perplexity, verbose = verbose)
    object <- Seurat::FindNeighbors(object, reduction = 'pca', k.param = k.param, dims = dims, verbose = verbose)
    object <- Seurat::FindClusters(object, algorithm = algorithm, resolution = resolution, verbose = verbose)
  }else if(assay=='RNA'){
    object <- Seurat::FindNeighbors(object, reduction = 'pca', k.param = k.param, dims = dims, verbose = verbose)
    object <- Seurat::FindClusters(object, algorithm = algorithm, resolution = resolution, verbose = verbose)
    object <- Seurat::RunUMAP(object, reduction = "pca", dims = dims, umap.method = 'umap-learn', metric = 'correlation', reduction.key = "UMAP_", reduction.name = "umap", verbose = verbose)
    object <- Seurat::RunTSNE(object, reduction = "pca", dims = dims, reduction.key = 'TSNE_', reduction.name = 'tsne', perplexity = perplexity, verbose = verbose)
  }
  return(object)
}
