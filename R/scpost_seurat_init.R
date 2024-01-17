#' Function to build a Seurat object, then normalize, scale & run PCA
#'
#' This function aggregates several Seurat package functions to build a Seurat object, then normalize data, scale data, & run PCA
#' @param counts Counts matrix
#' @param pc_max Max number of PCs to incorporate from PCA
#' @param verbose TRUE or FALSE; print Seurat package messages
#' @return Seurat object
#' @export
#' @examples
#' scpost_seurat_init(counts=counts, pc_max=10, verbose=TRUE)
#

#
scpost_seurat_init <- function(
    counts=counts,
    pc_max=10,
    verbose=TRUE) {

    # Initialize the Seurat object with the raw (non-normalized data)
    if (verbose == TRUE) { cat("   Seurat - Create", "\n") }
    seurat.obj <- Seurat::CreateSeuratObject(counts = counts,
                min.cells = 0,
                min.features = 0)
    #
    if (verbose == TRUE) { cat("   Seurat - Normalize", "\n") }
    seurat.obj <- Seurat::NormalizeData(seurat.obj,
                normalization.method = "LogNormalize",
                scale.factor = 10000,
                verbose = verbose)
    #    
    if (verbose == TRUE) { cat("   Seurat - Find Variable Features", "\n") }
    seurat.obj <- Seurat::FindVariableFeatures(seurat.obj,
                selection.method = "vst",
                nfeatures = 2000,
                verbose = verbose)
    #
    if (verbose == TRUE) { cat("   Seurat - Scale", "\n") }
    seurat.obj <- Seurat::ScaleData(seurat.obj,
                verbose=verbose)
    #
    if (verbose == TRUE) { cat("   Seurat - PCA", "\n") }
    seurat.obj <- Seurat::RunPCA(seurat.obj,
            features = Seurat::VariableFeatures(object = seurat.obj),
            verbose = verbose,
            npcs = pc_max,
            nfeatures.print = 2,
            ndims.print = 1:2)
    #
    return(seurat.obj)
    #
}













