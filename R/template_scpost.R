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
	cat("Read ExpressionSet", "\n")
	if (file.exists(file.path(dir_input, "ExpressionSet.rds"))) {
		dataset <- readRDS(file.path(dir_input, "ExpressionSet.rds"));
	} else {
		stop("No ExpressionSet Detected!")
	}

	#
	cat("Store Parameters in ExpressionSet assayData Params", "\n")
	Biobase::assayData(dataset)$Params["scpost_Parameters"] <- list(param.list)

	# Set 'reduction' to 'pca' to no batch correction
	if (param.list[["batch_correction_method"]] == FALSE) {
		reduction <- "pca"
	} else {
		reduction <- param.list[["batch_correction_method"]]
	}

	#
	cat("Select Cells", "\n")
	cells <- dataset$ID[which(dataset$Cell_Filter == "Cell")]
	#
	cat("Select Genes", "\n")
	genes <- Biobase::fData(dataset)$Ensembl[which(fData(dataset)$Gene_Filter == "Expressed")]

	#
	seurat.obj <- scpost::scpost_seurat_init(
    	counts=Biobase::exprs(dataset)[genes, cells],
    	pc_max=param.list[["pc_max"]],
    	verbose=verbose)
	#
	Biobase::assayData(dataset)$Seurat["PCA"] <- list(Seurat::Embeddings(seurat.obj, 'pca'))
	#
	saveRDS(Seurat::Embeddings(seurat.obj, 'pca'), file.path(file.path(dir_output, dir_analysis), "PCA.rds"))

	#
	if (param.list[["batch_var"]] != FALSE & reduction == "harmony") {
		#
		batch_metadata <- Biobase::pData(dataset)[colnames(seurat.obj), param.list[["batch_var"]]]
		names(batch_metadata) <- colnames(seurat.obj)
		#
		seurat.obj <- scpost::scpost_batch_correction(
    		seurat.obj=seurat.obj,
    		batch_metadata=batch_metadata,
    		method=reduction,
    		pc_max=param.list[["pc_max"]],
    		verbose=verbose)
		#
		Biobase::assayData(dataset)$Seurat["PCA_Harmony"] <- list(Seurat::Embeddings(seurat.obj, 'harmony'))
		#
		saveRDS(Seurat::Embeddings(seurat.obj, 'harmony'), file.path(file.path(dir_output, dir_analysis), "PCA_Harmony.rds"))
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
		Biobase::assayData(dataset)$Seurat[["Analysis"]][paste("Seurat_Res-", res, "_Clusters", sep="")] <- list(length(unique(clusters)))
	    #
	   	not.clusters <- rep(NA, length(setdiff(dataset$ID, names(clusters))))
	   	names(not.clusters) <- setdiff(dataset$ID, names(clusters))
	    seurat.clusters.all <- c(clusters, not.clusters)[dataset$ID]
	    #
	    cat(paste("Resolution", res, "- Clusters:", length(unique(clusters))), "\n")

	    # Store clusters in pData slot
		dataset[[paste("Seurat_Clusters_Res-", res, sep="")]] <- as.factor(seurat.clusters.all)
	}


    #
    pdata.out <- Biobase::pData(dataset)[,c("Sample", "Cell_Filter", grep("Seurat_Clusters_Res-", colnames(Biobase::pData(dataset)), value=TRUE))]
	#
	saveRDS(pdata.out, file.path(file.path(dir_output, dir_analysis), "MetaData_Seurat_Clusters.rds"))

	
	#
	if (param.list[["umap"]] == TRUE) {
		#
		dir_umap <- "UMAP"
		dir.create(file.path(file.path(dir_output, dir_analysis), dir_umap))
		#
		for (iter in 1:param.list[["iterations"]]) {
			#
    		cat(paste("Run UMAP - Iteration:", iter), "\n")
			set.seed(Biobase::assayData(dataset)$Params[["Seeds"]][iter])
        	seurat.obj <- suppressWarnings(Seurat::RunUMAP(
	    	    object=seurat.obj,
	    	    reduction=reduction,
	    	    dims=1:param.list[["pc_max"]],
	    	    seed.use=Biobase::assayData(dataset)$Params[["Seeds"]][iter],
	    	    verbose=verbose))

			# Store dimred output in ESet
	    	Biobase::assayData(dataset)$Seurat[["UMAP"]][paste("Iteration_", iter, sep="")] <- list(seurat.obj$umap@cell.embeddings)
			# Save dimred output as .rds file
			saveRDS(seurat.obj$umap@cell.embeddings, file.path(file.path(file.path(dir_output, dir_analysis), dir_umap), paste("UMAP_Iteration_", iter, ".rds", sep="")))
		}
	}

	#
	if (param.list[["tsne"]] == TRUE) {
		#
		dir_tsne <- "tSNE"
		dir.create(file.path(file.path(dir_output, dir_analysis), dir_tsne))
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
	    	    seed.use=Biobase::assayData(dataset)$Params[["Seeds"]][iter],
	    	   	nthreads=threads,
		        max_iter=2000,
		        check_duplicates=FALSE,
	    	    verbose=verbose))

			# Store dimred output in ESet
	    	Biobase::assayData(dataset)$Seurat[["tSNE"]][paste("Iteration_", iter, sep="")] <- list(seurat.obj$tsne@cell.embeddings)
			# Save dimred output as .rds file
			saveRDS(seurat.obj$tsne@cell.embeddings, file.path(file.path(file.path(dir_output, dir_analysis), dir_tsne), paste("tSNE_Iteration_", iter, ".rds", sep="")))
		}
	}

	# save ExpressionSet as RDS
	cat("Save ExpressionSet", "\n")
	saveRDS(dataset, file.path(dir_output, "ExpressionSet.rds"));
	#
	return(dataset)
	#
}




