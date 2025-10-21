#' Template function for post-processing scRNA-Seq data
#'
#' This function reads an ExpressionSet object with scRNA data and performs core functionalities for post-processing of scRNA-Seq data
#' @param dir_output Output directory for post-processing analysis results
#' @return ExpressionSet object
#' @export
#' @examples
#' template_scpost(dir_output="/path/to/output")
#

#
template_scpost <- function(
	dir_output=dir_output) {
	#
	cat(paste(Sys.time()), "\n")
	cat("Initiate scRNA-Seq Post-Postprocessing", "\n")

	# Read Parameters
	parameters <- read.csv(file.path(dir_output, "scpost_parameters.csv"), stringsAsFactors=FALSE, row.names=1)
	
	# Pass parameters to variables & incorporate into ExpressionSet
	param.list <- lapply(rownames(parameters), function(param) {
		#
		cat(paste(parameters[param, "Brief"], parameters[param, "Selection"]), "\n")
		return(parameters[param, "Selection"])
		#
	})
	#
	names(param.list) <- rownames(parameters)
	#
	param.list["pc_max"] <- list(as.numeric(param.list[["pc_max"]]))
	param.list["res_ops"] <- list(as.numeric(unlist(strsplit(param.list[["res_ops"]], ", "))))
	param.list["iterations"] <- list(as.numeric(param.list[["iterations"]]))
	param.list["q_max"] <- list(as.numeric(param.list[["q_max"]]))
	param.list["z_min"] <- list(as.numeric(param.list[["z_min"]]))
	
	#
	verbose <- TRUE
	dir_input <- param.list[["dir_input"]]
	dir_output <- param.list[["dir_output"]]
	batch_var <- param.list[["batch_var"]]
	reduction <- param.list[["batch_correction_method"]]

	#
	dir_analysis <- "scpost_analysis"
	dir.create(file.path(dir_output, dir_analysis))

	#
	cat("Read Dataset", "\n")
	# Try to detect object type by checking for different file names
	if (file.exists(file.path(dir_input, "ExpressionSet.rds"))) {
		dataset <- readRDS(file.path(dir_input, "ExpressionSet.rds"))
		object_type <- "eset"
		cat("Loaded ExpressionSet object", "\n")
	} else if (file.exists(file.path(dir_input, "seurat.rds"))) {
		dataset <- readRDS(file.path(dir_input, "seurat.rds"))
		object_type <- "seurat"
		cat("Loaded Seurat object", "\n")
	} else if (file.exists(file.path(dir_input, "sce.rds"))) {
		dataset <- readRDS(file.path(dir_input, "sce.rds"))
		object_type <- "sce"
		cat("Loaded SingleCellExperiment object", "\n")
	} else {
		stop("No dataset detected! Expected ExpressionSet.rds, seurat.rds, or sce.rds in input directory")
	}

	#
	cat("Store Parameters in dataset", "\n")
	if (object_type == "eset") {
		Biobase::assayData(dataset)$Params["scpost_Parameters"] <- list(param.list)
	} else if (object_type == "seurat") {
		dataset@misc$scpost_Parameters <- param.list
	} else if (object_type == "sce") {
		S4Vectors::metadata(dataset)$scpost_Parameters <- param.list
	}

	# Set 'reduction' to 'pca' to no batch correction
	if (param.list[["batch_correction_method"]] == FALSE) {
		reduction <- "pca"
	} else {
		reduction <- param.list[["batch_correction_method"]]
	}

	#
	cat("Select Cells and Genes", "\n")
	# Extract cells and genes based on object type
	if (object_type == "eset") {
		cells <- dataset$ID[which(dataset$Cell_Filter == "Cell")]
		genes <- Biobase::fData(dataset)$Ensembl[which(Biobase::fData(dataset)$Gene_Filter == "Expressed")]
	} else if (object_type == "seurat") {
		cells <- rownames(dataset@meta.data)[which(dataset@meta.data$Cell_Filter == "Cell")]
		# For Seurat, use all genes since gene filtering is typically handled differently
		genes <- rownames(dataset)
		cat("Note: Using all genes from Seurat object (gene filtering not applied)", "\n")
	} else if (object_type == "sce") {
		cells <- rownames(SummarizedExperiment::colData(dataset))[which(SummarizedExperiment::colData(dataset)$Cell_Filter == "Cell")]
		# For SCE, use all genes since gene filtering is typically handled differently
		genes <- rownames(dataset)
		cat("Note: Using all genes from SingleCellExperiment object (gene filtering not applied)", "\n")
	}

	#
	# Extract counts matrix based on object type
	if (object_type == "eset") {
		counts_matrix <- Biobase::exprs(dataset)[genes, cells]
	} else if (object_type == "seurat") {
		# Subset Seurat object to selected cells only
		seurat.obj <- subset(dataset, cells = cells)
		counts_matrix <- NULL  # Will pass seurat object directly
	} else if (object_type == "sce") {
		counts_matrix <- SummarizedExperiment::assay(dataset, "counts")[genes, cells]
	}

	# Initialize Seurat processing
	if (object_type == "seurat") {
		# Seurat object already exists, just need to process it
		cat("Processing existing Seurat object", "\n")
		# Run normalization and PCA steps
		if (verbose == TRUE) { cat("   Seurat - Normalize", "\n") }
		seurat.obj <- Seurat::NormalizeData(seurat.obj,
			normalization.method = "LogNormalize",
			scale.factor = 10000,
			verbose = verbose)
		if (verbose == TRUE) { cat("   Seurat - Find Variable Features", "\n") }
		seurat.obj <- Seurat::FindVariableFeatures(seurat.obj,
			selection.method = "vst",
			nfeatures = 2000,
			verbose = verbose)
		if (verbose == TRUE) { cat("   Seurat - Scale", "\n") }
		seurat.obj <- Seurat::ScaleData(seurat.obj, verbose=verbose)
		if (verbose == TRUE) { cat("   Seurat - PCA", "\n") }
		seurat.obj <- Seurat::RunPCA(seurat.obj,
			features = Seurat::VariableFeatures(object = seurat.obj),
			verbose = verbose,
			npcs = param.list[["pc_max"]],
			nfeatures.print = 2,
			ndims.print = 1:2)
	} else {
		# For ExpressionSet and SCE, use scpost_seurat_init
		seurat.obj <- scpost::scpost_seurat_init(
			counts=counts_matrix,
			pc_max=param.list[["pc_max"]],
			verbose=verbose)
	}
	#
	# Store PCA results based on object type
	pca_embeddings <- Seurat::Embeddings(seurat.obj, 'pca')
	if (object_type == "eset") {
		Biobase::assayData(dataset)$Seurat["PCA"] <- list(pca_embeddings)
	} else if (object_type == "seurat") {
		# PCA already stored in Seurat object reductions slot
		cat("PCA stored in Seurat object @reductions$pca", "\n")
	} else if (object_type == "sce") {
		SingleCellExperiment::reducedDim(dataset, "PCA") <- pca_embeddings
	}
	#
	saveRDS(pca_embeddings, file.path(file.path(dir_output, dir_analysis), "PCA.rds"))

	#
	if (param.list[["batch_var"]] != FALSE & reduction == "harmony") {
		#
		# Extract batch metadata based on object type
		if (object_type == "eset") {
			batch_metadata <- Biobase::pData(dataset)[colnames(seurat.obj), param.list[["batch_var"]]]
		} else if (object_type == "seurat") {
			batch_metadata <- dataset@meta.data[colnames(seurat.obj), param.list[["batch_var"]]]
		} else if (object_type == "sce") {
			batch_metadata <- SummarizedExperiment::colData(dataset)[colnames(seurat.obj), param.list[["batch_var"]]]
		}
		names(batch_metadata) <- colnames(seurat.obj)
		#
		seurat.obj <- scpost::scpost_batch_correction(
    		seurat.obj=seurat.obj,
    		batch_metadata=batch_metadata,
    		method=reduction,
    		pc_max=param.list[["pc_max"]],
    		verbose=verbose)
		#
		# Store harmony results based on object type
		harmony_embeddings <- Seurat::Embeddings(seurat.obj, 'harmony')
		if (object_type == "eset") {
			Biobase::assayData(dataset)$Seurat["PCA_Harmony"] <- list(harmony_embeddings)
		} else if (object_type == "seurat") {
			# Harmony already stored in Seurat object reductions slot
			cat("Harmony stored in Seurat object @reductions$harmony", "\n")
		} else if (object_type == "sce") {
			SingleCellExperiment::reducedDim(dataset, "PCA_Harmony") <- harmony_embeddings
		}
		#
		saveRDS(harmony_embeddings, file.path(file.path(dir_output, dir_analysis), "PCA_Harmony.rds"))
	}

	#
	for (res in param.list[["res_ops"]]) {
		#
		cat(paste("Run Seurat Clustering - Resolution:", res), "\n")
		#
		seurat.obj <- scpost::scpost_seurat_cluster(
	    	seurat.obj=seurat.obj,
	    	pc_max=param.list[["pc_max"]],
	    	resolution=res,
	    	reduction=reduction,
	    	seurat_return=TRUE,
	    	verbose=verbose)
	    #
	    clusters <- as.numeric(Seurat::Idents(seurat.obj))
	    names(clusters) <- colnames(seurat.obj)

	    # Store numbers of clusters per resolution value
	    if (object_type == "eset") {
			Biobase::assayData(dataset)$Seurat[["Analysis"]][paste("Seurat_Res-", res, "_Clusters", sep="")] <- list(length(unique(clusters)))
		} else if (object_type == "seurat") {
			if (is.null(dataset@misc$Seurat_Analysis)) {
				dataset@misc$Seurat_Analysis <- list()
			}
			dataset@misc$Seurat_Analysis[[paste("Seurat_Res-", res, "_Clusters", sep="")]] <- length(unique(clusters))
		} else if (object_type == "sce") {
			if (is.null(S4Vectors::metadata(dataset)$Seurat_Analysis)) {
				S4Vectors::metadata(dataset)$Seurat_Analysis <- list()
			}
			S4Vectors::metadata(dataset)$Seurat_Analysis[[paste("Seurat_Res-", res, "_Clusters", sep="")]] <- length(unique(clusters))
		}
	    #
	    # Get all cell IDs based on object type
	    if (object_type == "eset") {
	    	all_cell_ids <- dataset$ID
	    } else if (object_type == "seurat") {
	    	all_cell_ids <- rownames(dataset@meta.data)
	    } else if (object_type == "sce") {
	    	all_cell_ids <- rownames(SummarizedExperiment::colData(dataset))
	    }

	   	not.clusters <- rep(NA, length(setdiff(all_cell_ids, names(clusters))))
	   	names(not.clusters) <- setdiff(all_cell_ids, names(clusters))
	    seurat.clusters.all <- c(clusters, not.clusters)[all_cell_ids]
	    #
	    cat(paste("Resolution", res, "- Clusters:", length(unique(clusters))), "\n")

	    # Store clusters in metadata based on object type
	    if (object_type == "eset") {
			dataset[[paste("Seurat_Clusters_Res-", res, sep="")]] <- as.factor(seurat.clusters.all)
		} else if (object_type == "seurat") {
			dataset@meta.data[[paste("Seurat_Clusters_Res-", res, sep="")]] <- as.factor(seurat.clusters.all)
		} else if (object_type == "sce") {
			SummarizedExperiment::colData(dataset)[[paste("Seurat_Clusters_Res-", res, sep="")]] <- as.factor(seurat.clusters.all)
		}
	}


    #
    # Extract metadata for output based on object type
    if (object_type == "eset") {
    	pdata.out <- Biobase::pData(dataset)[,c("Sample", "Cell_Filter", grep("Seurat_Clusters_Res-", colnames(Biobase::pData(dataset)), value=TRUE))]
    } else if (object_type == "seurat") {
    	cluster_cols <- grep("Seurat_Clusters_Res-", colnames(dataset@meta.data), value=TRUE)
    	pdata.out <- dataset@meta.data[,c("Sample", "Cell_Filter", cluster_cols)]
    } else if (object_type == "sce") {
    	cluster_cols <- grep("Seurat_Clusters_Res-", colnames(SummarizedExperiment::colData(dataset)), value=TRUE)
    	pdata.out <- as.data.frame(SummarizedExperiment::colData(dataset)[,c("Sample", "Cell_Filter", cluster_cols)])
    }
	#
	saveRDS(pdata.out, file.path(file.path(dir_output, dir_analysis), "MetaData_Seurat_Clusters.rds"))

	
	#
	if (param.list[["umap"]] == TRUE) {
		#
		dir_umap <- "UMAP"
		dir.create(file.path(file.path(dir_output, dir_analysis), dir_umap))
		#
		# Get seeds based on object type
		if (object_type == "eset") {
			seeds <- Biobase::assayData(dataset)$Params[["Seeds"]]
		} else if (object_type == "seurat") {
			# Check if seeds exist in misc slot, otherwise generate new ones
			if (is.null(dataset@misc$Seeds)) {
				dataset@misc$Seeds <- scprep::scprep_seeds(n.seeds=1000)
			}
			seeds <- dataset@misc$Seeds
		} else if (object_type == "sce") {
			# Check if seeds exist in metadata, otherwise generate new ones
			if (is.null(S4Vectors::metadata(dataset)$Seeds)) {
				S4Vectors::metadata(dataset)$Seeds <- scprep::scprep_seeds(n.seeds=1000)
			}
			seeds <- S4Vectors::metadata(dataset)$Seeds
		}
		#
		for (iter in 1:param.list[["iterations"]]) {
			#
    		cat(paste("Run UMAP - Iteration:", iter), "\n")
			set.seed(seeds[iter])
        	seurat.obj <- suppressWarnings(Seurat::RunUMAP(
	    	    object=seurat.obj,
	    	    reduction=reduction,
	    	    dims=1:param.list[["pc_max"]],
	    	    seed.use=seeds[iter],
	    	    verbose=verbose))

			# Store dimred output based on object type
			umap_embeddings <- seurat.obj$umap@cell.embeddings
			if (object_type == "eset") {
	    		Biobase::assayData(dataset)$Seurat[["UMAP"]][paste("Iteration_", iter, sep="")] <- list(umap_embeddings)
	    	} else if (object_type == "seurat") {
	    		# Store in misc slot for non-active UMAP iterations
	    		if (is.null(dataset@misc$UMAP_iterations)) {
	    			dataset@misc$UMAP_iterations <- list()
	    		}
	    		dataset@misc$UMAP_iterations[[paste("Iteration_", iter, sep="")]] <- umap_embeddings
	    	} else if (object_type == "sce") {
	    		# Store in reducedDims with iteration-specific names
	    		SingleCellExperiment::reducedDim(dataset, paste("UMAP_Iteration_", iter, sep="")) <- umap_embeddings
	    	}
			# Save dimred output as .rds file
			saveRDS(umap_embeddings, file.path(file.path(file.path(dir_output, dir_analysis), dir_umap), paste("UMAP_Iteration_", iter, ".rds", sep="")))
		}
	}

	#
	if (param.list[["tsne"]] == TRUE) {
		#
		dir_tsne <- "tSNE"
		dir.create(file.path(file.path(dir_output, dir_analysis), dir_tsne))
		#
		# Get seeds based on object type (same as UMAP)
		if (object_type == "eset") {
			seeds <- Biobase::assayData(dataset)$Params[["Seeds"]]
		} else if (object_type == "seurat") {
			# Check if seeds exist in misc slot, otherwise generate new ones
			if (is.null(dataset@misc$Seeds)) {
				dataset@misc$Seeds <- scprep::scprep_seeds(n.seeds=1000)
			}
			seeds <- dataset@misc$Seeds
		} else if (object_type == "sce") {
			# Check if seeds exist in metadata, otherwise generate new ones
			if (is.null(S4Vectors::metadata(dataset)$Seeds)) {
				S4Vectors::metadata(dataset)$Seeds <- scprep::scprep_seeds(n.seeds=1000)
			}
			seeds <- S4Vectors::metadata(dataset)$Seeds
		}
		#
		for (iter in 1:param.list[["iterations"]]) {
			#
    		cat(paste("Run tSNE - Iteration:", iter), "\n")
			threads <- parallel::detectCores()
	    	cat(paste("Available Cores:", threads), "\n")
        	seurat.obj <- suppressWarnings(Seurat::RunTSNE(
	    	    object=seurat.obj,
	    	    reduction=reduction,
	    	    dims=1:param.list[["pc_max"]],
	    	    seed.use=seeds[iter],
	    	   	nthreads=threads,
		        max_iter=2000,
		        check_duplicates=FALSE,
	    	    verbose=verbose))

			# Store dimred output based on object type
			tsne_embeddings <- seurat.obj$tsne@cell.embeddings
			if (object_type == "eset") {
	    		Biobase::assayData(dataset)$Seurat[["tSNE"]][paste("Iteration_", iter, sep="")] <- list(tsne_embeddings)
	    	} else if (object_type == "seurat") {
	    		# Store in misc slot for non-active tSNE iterations
	    		if (is.null(dataset@misc$tSNE_iterations)) {
	    			dataset@misc$tSNE_iterations <- list()
	    		}
	    		dataset@misc$tSNE_iterations[[paste("Iteration_", iter, sep="")]] <- tsne_embeddings
	    	} else if (object_type == "sce") {
	    		# Store in reducedDims with iteration-specific names
	    		SingleCellExperiment::reducedDim(dataset, paste("tSNE_Iteration_", iter, sep="")) <- tsne_embeddings
	    	}
			# Save dimred output as .rds file
			saveRDS(tsne_embeddings, file.path(file.path(file.path(dir_output, dir_analysis), dir_tsne), paste("tSNE_Iteration_", iter, ".rds", sep="")))
		}
	}

	# Save dataset based on object type
	if (object_type == "eset") {
		cat("Save ExpressionSet", "\n")
		saveRDS(dataset, file.path(dir_output, "ExpressionSet.rds"))
	} else if (object_type == "seurat") {
		cat("Save Seurat object", "\n")
		saveRDS(dataset, file.path(dir_output, "seurat.rds"))
	} else if (object_type == "sce") {
		cat("Save SingleCellExperiment object", "\n")
		saveRDS(dataset, file.path(dir_output, "sce.rds"))
	}
	#
	return(dataset)
	#
}




