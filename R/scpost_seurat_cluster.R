#' Function to perform Seurat clustering using a Seurat object generated with scpost_seurat_init
#'
#' This function applies Seurat package functions to a Seurat object in order to find nearest neighbors & perform SNN & louvain community detection
#' @param seurat.obj A Seurat object generated with scpost_seurat_init
#' @param pc_max Number of PCs to utilize for clustering
#' @param resolution Numeric "resolution" value that dictates cluster size
#' @param reduction Default is "pca"; Select "harmony" if a batch_correction variable was specified on scpost_parameters.csv (batch_correction != FALSE)
#' @param seurat_return If TRUE, returns a Seurat object; if FALSE, returns cluster assignments
#' @param verbose TRUE or FALSE; print Seurat package messages
#' @return Seurat object or cluster assignments
#' @export
#' @examples
#' scpost_seurat_cluster(seurat.obj=seurat.obj, npc=10, resolution=0.05, reduction="pca", seurat_return=FALSE, verbose=TRUE)
#

#
scpost_seurat_cluster <- function(
    seurat.obj=seurat.obj,
    pc_max=10,
    resolution=0.05,
    reduction="pca",
    seurat_return=FALSE,
    verbose=TRUE) {
    #
    if (verbose == TRUE) { cat("   Seurat - Nearest Neighbors", "\n") }
    seurat.obj <- Seurat::FindNeighbors(object = seurat.obj,
            reduction = reduction,
            dims = 1:pc_max,
            verbose=verbose)
    #
    if (verbose == TRUE) { cat("   Seurat - Shared Nearest Neighbors & Louvain Community Detection", "\n") }
    seurat.obj <- Seurat::FindClusters(seurat.obj,
        resolution = resolution,
        algorithm = 1,
        verbose=verbose,
        n.start = 10)
    #
    if (verbose == TRUE) { cat("   Seurat - Complete", "\n") }
    #
    if (seurat_return == TRUE) {
        return(seurat.obj)
    } else {
        return(Idents(seurat.obj))
    }
}
