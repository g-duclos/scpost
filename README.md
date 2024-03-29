![Logo](inst/extdata/scpost_Logo.png)

***

# scpost
An R package for post-processing single-cell RNA-Seq data involving unsupervised clustering, dimensionality reduction, batch correction, sample-cluster composition calculation, pseudobulk differential expression (in progress), and functional annotation (in progress). This package is designed to be used to process data that has been aggregated by the ["scprep" R package](https://github.com/g-duclos/scprep).

***

## Installation

The **scpost** R package can be installed from Github using devtools:
```
devtools::install_github("g-duclos/scpost")
```

***

## Getting Started

#### Parameters
Specify the pipeline parameters - view the parameters file here: [scpost_parameters.csv](inst/extdata/scpost_parameters.csv) (click "Raw file" at the top right to download)

**Critical Parameters**
* *dir_input* (path to input directory)
* *dir_output* (path to output directory)

#### Input Directory
Define the input directory (*dir_input*), which must contain the following file generated by the ["scprep" R package](https://github.com/g-duclos/scprep):
* **ExpressionSet.rds** (S4 object containing GEX counts + sample & feature metadata)

#### Dependency Note:
The ["Seurat" R package](https://satijalab.org/seurat/) must be installed (v3, v4, or v5 is acceptable) in order to use **scpost**. However, Seurat is NOT included as a formal package dependency due to common installation complications.

***

## Overview

Template function to perform unsupervised clustering, dimensionality reduction, batch correction, sample-cluster composition calculation, pseudobulk differential expression (in progress), and functional annotation (in progress).

```
library(Biobase)

dataset <- scpost::template_scpost(dir_output=dir_output)
```

**Core functions:**

<details>
	<summary>
		A single command to add GEX counts data to a Seurat object, then perform standard normalization (natural log counts per ten thousand), select the top 2000 highly variable genes, scale and center the normalized data, and perform principal component analysis (PCA).
	</summary>
<pre>
# Select high quality cells
cells <- dataset$ID[which(dataset$Cell_Filter == "Cell")]
</pre>
<pre>
# Select genes with minimum expression requirements
genes <- Biobase::fData(dataset)$Ensembl[which(fData(dataset)$Gene_Filter == "Expressed")]
</pre>
<pre>
# Build Seurat object, normalize, select HVG, scale, & PCA
seurat.obj <- scpost::scpost_seurat_init(
    counts=counts[genes, cells],
    pc_max=pc_max,
    verbose=verbose)
</pre>
</details>


<details>
	<summary>
		Data integration intended to correct for technical batch effects for downstream clustering & dimensionality reduction, based on a selected "batch" variable of interest using one of the following methods:
		<ul><li>
			<a href="https://github.com/immunogenomics/harmony">Harmony R package</a>
		</li>
		<li>
			<a href="https://satijalab.org/seurat/articles/integration_rpca.html">Seurat's RPCA method</a>
		</li>
	</summary>
<pre>
# Select batch metadata according to metadata variable of interest
batch_metadata <- Biobase::pData(dataset)[colnames(seurat.obj), batch_var]
names(batch_metadata) <- colnames(seurat.obj)
</pre>
<pre>
# Perform batch correction with "harmony" or "rpca" method
seurat.obj <- scpost::scpost_batch_correction(
   	seurat.obj=seurat.obj,
    batch_metadata=batch_metadata,
    method=reduction,
    pc_max=pc_max,
    verbose=verbose)
</pre>
</details>

<details>
	<summary>
		A single command to perform unsupervised cell clustering using Seurat's standard methodology (SNN & Louvain community detection) with the option of returning a range of cluster solutions.
	</summary>
<pre>
# Perform unsupervised clustering
seurat.obj <- scpost::scpost_seurat_cluster(
	seurat.obj=seurat.obj,
	pc_max=pc_max,
	resolution=res,
	reduction=reduction,
	seurat_return=TRUE,
	verbose=verbose)
</pre>
</details>


<details>
	<summary>
	The Seurat *RunTSNE* function is used to run t-SNE on the dataset.
	</summary>
	<ul><li>
		tSNE results are stored as a list labeled "tSNE" in the "assayData" slot labeled "Seurat" in the ExpressionSet object.
	<li>
		tSNE results will reflect the decision to use or not use one of the specified batch correction methods
	</li>
	<li>
		According to the number of iterations specified in the "scpost_parameters.csv" file, results are subsequently stored as "Iteration_1", "Iteration_2", ...
	</li>
	<li>
		According to the number of iterations specified in the "scpost_parameters.csv" file, results are saved in a subdirectory named "scpost_analysis" in *dir_output* in the RDS format as "tSNE_Iteration_1.rds", "tSNE_Iteration_2.rds", ...
	</li>
	</ul>
</details>


<details>
	<summary>
		The Seurat *RunUMAP* function is used to run UMAP on the dataset and results/
	</summary>
	<ul><li>
		UMAP results are stored as a list labeled "UMAP" in the "assayData" slot labeled "Seurat" in the ExpressionSet object.
	</li>
	<li>
		UMAP results will reflect the decision to use or not use one of the specified batch correction methods
	</li>
	<li>
		According to the number of iterations specified in the "scpost_parameters.csv" file, results are subsequently stored as "Iteration_1", "Iteration_2", ...
	</li>
	<li>
		According to the number of iterations specified in the "scpost_parameters.csv" file, results are saved in a subdirectory named "scpost_analysis" in *dir_output* in the RDS format as "UMAP_Iteration_1.rds", "UMAP_Iteration_2.rds", ...
	</li>
	</ul>
</details>


<li>
	Sample-cluster composition assessment methods with built-in statistical methods for univariate or multivariate analysis of cluster-associations with variable(s) of interest are a work in progress ...
</li>
<br>


<li>
	Pseudobulk differential expression functionalities, involving the use of the <a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">DESeq2 R package</a> are a work in progress...
</li>
<br>


<li>
	Functional annotation, involving the use of gene set databases, including <a href="https://www.gsea-msigdb.org/gsea/msigdb/">the Molecular Signatures Database (MSigDB)</a> and <a href="https://geneontology.org">Gene Ontology (GO</a>, with methods, including the <a href="https://bioconductor.org/packages/release/bioc/html/fgsea.html">Fast Gene Set Enrichment analysis (fgsea) R package</a> and the <a href="https://cran.r-project.org/web/packages/enrichR/index.html">enrichR R package</a>
</li>
<br>


<details>
	<summary>Additional features</summary>
<ul><li>
	Store parameters specified in the scpost_parameters.csv file in a list labeled "Parameters" in the "assayData" slot labeled "Params" in the ExpressionSet object
</li>

<li>
	Save ExpressionSet RDS file in *dir_output*
</li>
</ul>
</details>

***
