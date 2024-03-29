% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scpost_batch_correction.R
\name{scpost_batch_correction}
\alias{scpost_batch_correction}
\title{Function to perform batch correction for scRNA-Seq data}
\usage{
scpost_batch_correction(
  seurat.obj = seurat.obj,
  batch_metadata = batch_metadata,
  method = "harmony",
  pc_max = 10,
  verbose = TRUE
)
}
\arguments{
\item{seurat.obj}{A Seurat object generated with scpost_seurat_init}

\item{batch_metadata}{Vector containing categorical batch specification for each cell; vector length MUST be equivalent to number of cells in the Seurat object; vector elements MUST be labeled/named with cell IDs}

\item{method}{Options for batch correction methods, including "harmony" or "cca" (original Seurat method) (both are compatible with Seurat objects)}

\item{pc_max}{Max number of PCs to incorporate for PCA - this should be equivalent to the number of PCs you will use as input for downstream clustering and dimensionality reduction}

\item{verbose}{TRUE or FALSE; print package messages}
}
\value{
Seurat object containing batch correction results (in the "reduction" slot for Harmony, assay set to "integrated" for CCA)
}
\description{
This function  reads in a Seurat object and corrects for batch effects according to a specified 'batch' variable
}
\examples{
scpost_batch_correction(seurat.obj=seurat.obj, batch_metadata=batch_metadata, method="harmony", pc_max=10, verbose=TRUE)
}
