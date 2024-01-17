#' Function to perform batch correction for scRNA-Seq data
#'
#' This function  reads in a Seurat object and corrects for batch effects according to a specified 'batch' variable
#' @param seurat.obj A Seurat object generated with scpost_seurat_init
#' @param batch_metadata Vector containing categorical batch specification for each cell; vector length MUST be equivalent to number of cells in the Seurat object; vector elements MUST be labeled/named with cell IDs
#' @param method Options for batch correction methods, including "harmony" or "cca" (original Seurat method) (both are compatible with Seurat objects)
#' @param pc_max Max number of PCs to incorporate for PCA - this should be equivalent to the number of PCs you will use as input for downstream clustering and dimensionality reduction
#' @param verbose TRUE or FALSE; print package messages
#' @return Seurat object containing batch correction results (in the "reduction" slot for Harmony, assay set to "integrated" for CCA)
#' @export
#' @examples
#' scpost_batch_correction(seurat.obj=seurat.obj, batch_metadata=batch_metadata, method="harmony", pc_max=10, verbose=TRUE)
#

#
scpost_batch_correction <- function(
    seurat.obj=seurat.obj,
    batch_metadata=batch_metadata,
    method="harmony",
    pc_max=10,
    verbose=TRUE) {  
	#
    if (length(intersect(names(batch_metadata), rownames(seurat.obj@meta.data))) == nrow(seurat.obj@meta.data)) {
        #
        cat("Add Batch MetaData to Seurat Object", "\n")
        #
        seurat.obj@meta.data[["batch"]] <- batch_metadata[rownames(seurat.obj@meta.data)]
        #
    } else {
        cat("Batch MetaData does NOT match Seurat Object !!!", "\n")
        cat("Check batch_metadata vector length & names", "\n")
    }
  	#
  	if (method == "harmony") {
    	#
    	cat("Run 'Harmony' Batch Correction", "\n")
    	#
    	seurat.obj <- suppressWarnings(harmony::RunHarmony(seurat.obj, "batch",
    		dims.use = 1:pc_max,
    		verbose = verbose))
    	#
  	} else if (method == "cca") {
    	#
    	cat("Run 'RPCA' Batch Correction", "\n")
    	# Split seurat object by batch
    	seurat.obj.list <- Seurat::SplitObject(seurat.obj, split.by = "batch")
    	# Normalize & Find HVG
    	seurat.obj.list <- lapply(1:length(seurat.obj.list), function(i) {
    		#
    		cat("Batch", i, "of", length(seurat.obj.list), "-", names(seurat.obj.list)[i], "\n")
    		#
		    if (verbose == TRUE) { cat("   Seurat - Normalize", "\n") }
		    obj <- Seurat::NormalizeData(seurat.obj.list[[i]],
                normalization.method = "LogNormalize",
                scale.factor = 10000,
                verbose = verbose)
    		#    
    		if (verbose == TRUE) { cat("   Seurat - Scale", "\n") }
    		obj <- Seurat::FindVariableFeatures(obj,
                selection.method = "vst",
                nfeatures = 2000,
                verbose = verbose)
    		#
    		return(obj)
    		#
    	})
		# Select features for integration that are variable across all datasets
		if (verbose == TRUE) { cat("   Seurat - Select Integration Features", "\n") }
		features <- Seurat::SelectIntegrationFeatures(object.list = seurat.obj.list,
			verbose = verbose)
		# Run PCA on each dataset based on selected features
    	seurat.obj.list <- lapply(1:length(seurat.obj.list), function(i) {
	 		#
    		cat("Batch", i, "of", length(seurat.obj.list), "\n")
			#
		    if (verbose == TRUE) { cat("   Seurat - Scale", "\n") }
    		obj <- Seurat::ScaleData(seurat.obj.list[[i]],
    			features = features,
    			verbose = verbose)
    		#
    		if (verbose == TRUE) { cat("   Seurat - PCA", "\n") }
   			obj <- Seurat::RunPCA(obj,
    			features = features,
    			npcs = pc_max,
        		nfeatures.print = 2,
        		ndims.print = 1:2,
        		verbose = verbose)
    		#
    		return(obj)
    		#
		})
		#
		if (verbose == TRUE) { cat("   Seurat - Select Integration Features", "\n") }
		anchors <- Seurat::FindIntegrationAnchors(object.list = seurat.obj.list,
			anchor.features = features,
			reduction = "cca",
			verbose = verbose)
		# IntegrateData() creates an 'integrated' assay slot
		if (verbose == TRUE) { cat("   Seurat - Integrate Data", "\n") }
		seurat.obj <- Seurat::IntegrateData(anchorset = anchors,
			verbose = verbose)

		# Specify that downstream analysis is performed on the RPCA corrected data
		if (verbose == TRUE) { cat("   Seurat - Set Default Assay to 'integrated'", "\n") }
		Seurat::DefaultAssay(seurat.obj) <- "integrated"

		# Run the standard workflow for visualization and clustering
		if (verbose == TRUE) { cat("   Seurat - Scale Integrated Data", "\n") }
		seurat.obj <- Seurat::ScaleData(seurat.obj,
			verbose = verbose)
		#
		if (verbose == TRUE) { cat("   Seurat - Run PCA on Integrated Data", "\n") }
		seurat.obj <- Seurat::RunPCA(seurat.obj,
    		npcs = pc_max,
        	nfeatures.print = 2,
        	ndims.print = 1:2,
        	verbose = verbose)
		#
	  	return(seurat.obj)
		#	
  	} else {
    	#
    	stop("Accepted method NOT chosen - check spelling !!!")
    	#
  	}
}
