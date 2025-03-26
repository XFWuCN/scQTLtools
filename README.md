# scQTLtools:An R package for single-cell eQTL analysis.


## Introduction
Expression quantitative trait loci (eQTL) analysis links variations in 
gene expression levels to genotypes. This package attempts to identify genetic
variants that affect the expression of genes at a single-cell level, and can 
also do cis-eQTL analysis, visualize the results. 

## Citation

If you find this tool useful, please cite:

------------------------------------------------------------------------

***[https://github.com/XFWuCN/scQTLtools](https://github.com/XFWuCN/scQTLtools)***

***[https://bioconductor.org/packages/3.21/bioc/html/scQTLtools.html](https://bioconductor.org/packages/3.21/bioc/html/scQTLtools.html)***


------------------------------------------------------------------------

## Installation

```{r, eval = FALSE}
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("scQTLtools")
```


## Overview of the package

The functions in scQTLtools can be categorized into data input, data 
pre-process, sc-eQTL calling and visualization modules. Each of these functions
and a short description is summarized as shown below.

### General Workflow

Each module is summarized as shown below.

***[Overview](vignettes/Overview.jpg)***

scQTLtools requires two key input matrices: a genotype matrix and a gene 
expression matrix. The gene expression data object can be provided as a Seurat
v4 or a Bioconductor SingleCellExperiment object. It enables the analysis of genotype matrices in two forms: (1) ref and alt, and (2) ref/ref, alt/alt, and ref/alt. Additionally, the package includes functionality to filter Gene-SNP
pairs that are closely positioned in the genome, as nearby SNPs are more likely 
to influence gene expression. Moreover, visualization at the single-cell level demonstrates the specificity of eQTLs across distinct cell types or cellular 
states.

## Comparison and advantages compared to similar works

We compared scQTLtools to other packages with similar functionality, including 
eQTLsingle, SCeQTL, MatrixEQTL, and iBMQ, as shown in the table below.

***[Comparison](vignettes/Comparison.svg)***

Among these tools, scQTLtools stands out for its comprehensive features:

(1) scQTLtools accepts SingleCellExperiment objects and Seurat objects as
input data formats, which are particularly beneficial for users working with
single-cell RNA-seq data, and promote the interoperability with the current
Bioconductor ecosystem. 

(2) scQTLtools supports both binary and triple classification genotype 
matrices, enhancing its applicability across different genetic studies. 

(3) scQTLtools offers extensive data pre-processing capabilities, including 
quality control filtering for SNPs and genes, normalization of expression data, 
and customization of SNP-gene pair distances. This ensures high-quality and 
well-prepared input data for subsequent analysis. 

(4) scQTLtools provides three kinds of fitting models, which cater to various 
data distributions and analysis needs. This diversity allows users to select 
the most appropriate model for their specific dataset. 

(5) scQTLtools includes a range of visualization tools, these options 
facilitate detailed exploration and interpretation of eQTL results at the 
single-cell level. 

(6) scQTLtools supports grouping by both cell type and cell state, which is 
crucial for analyzing the nuanced effects of genetic variants on gene 
expression within heterogeneous cell populations.

Overall, scQTLtools offers a comprehensive suite of features that enhance the 
analysis and interpretation of eQTLs.

## Required input files
The input file requires genotype data, as well as a gene expression matrix or 
a Seurat object, or a SingleCellExperiment object.

- gene expression matrix: describes gene expressions, the row names represent
gene IDs or SYMBOL and the column names represent cell IDs.
- Seurat object: a Seurat object.
- SingleCellExperiment object: a SingleCellExperiment object.
- genotype matrix: A genotype matrix where each row is one variant and each 
column is one sample, and the scoring method is 0/1/2/3, 0 represents missing
values, 1 represents ref/ref, 2 represents alt/alt, and 3 represents ref/alt.

The columns of the genotype matrix should correspond to the columns of the gene
expression matrix.

**Example**

```{r input, message=FALSE}
library(scQTLtools)
# gene expression matrix
data(testGene)
# SeuratObject
data(testSeurat)
# load the genotype data
data(testSNP)
data(testSNP2)
```


## Create eqtl object
The createQTLObject class is an R object designed to store data related to eQTL
analysis, encompassing data lists, result data frames, and slots for 
biClassify, species, and group information.

**Example**

```{r createObject_matrix, message=FALSE}
eqtl_matrix <- createQTLObject(
    snpMatrix = testSNP,
    genedata = testGene,
    biClassify = FALSE,
    species = 'human',
    group = NULL)
```

Users can set biClassify to TRUE to change the genotype coding method.

**Example**

```{r createObject_matrix_bi, message=FALSE}
eqtl_matrix_bi <- createQTLObject(
    snpMatrix = testSNP,
    genedata = testGene,
    biClassify = TRUE,
    species = 'human',
    group = NULL)
```

Users can use Seuratobjct instead of gene expression matrix. 

**Example**

```{r createObject_seuratobject, message=FALSE}
eqtl_seurat <- createQTLObject(
    snpMatrix = testSNP2,
    genedata = testSeurat,
    biClassify = FALSE,
    species = 'human',
    group = "celltype")
```


## Normalize gene expression matrix
Use `normalizeGene()` to normalize the gene expression matrix.

**Example**

```{r Normalize_matrix, message=FALSE}
eqtl_matrix  <- normalizeGene(
    eQTLObject = eqtl_matrix, 
    method = "logNormalize")
```


## Identify the valid gene snp pairs
Here we use `filterGeneSNP()` to filter snp gene pairs.

**Example**

```{r filter_matrix, message=FALSE}
eqtl_matrix <- filterGeneSNP(
    eQTLObject = eqtl_matrix,
    snpNumOfCellsPercent = 2,
    expressionMin = 0,
    expressionNumOfCellsPercent = 2)
```

```{r filter_seuratobject, message=FALSE}
eqtl_seurat <- filterGeneSNP(
    eQTLObject = eqtl_seurat,
    snpNumOfCellsPercent = 2,
    expressionMin = 0,
    expressionNumOfCellsPercent = 2)
```


## Call single cell eQTL
Here we use `callQTL()` to do single cell eQTL analysis.

**Example**

```{r callQTL1_matrix, message=FALSE}
eqtl1_matrix <- callQTL(
    eQTLObject = eqtl_matrix,
    gene_ids = NULL,
    downstream = NULL,
    upstream = NULL,
    pAdjustMethod = "bonferroni",
    useModel = "poisson",
    pAdjustThreshold = 0.05,
    logfcThreshold = 0.1)
```

```{r callQTL1_seuratobject, message=FALSE}
eqtl1_seurat <- callQTL(
    eQTLObject = eqtl_seurat,
    gene_ids = NULL,
    downstream = NULL,
    upstream = NULL,
    pAdjustMethod = "bonferroni",
    useModel = "linear",
    pAdjustThreshold = 0.05,
    logfcThreshold = 0.025)
```

Users can use the parameter `gene_ids` to select one or several genes of 
interest for identifying sc-eQTLs.

**Example**

```{r callQTL2_matrix, message=FALSE}
eqtl2_matrix <- callQTL(
    eQTLObject = eqtl_matrix,
    gene_ids = c("CNN2", 
                "RNF113A", 
                "SH3GL1", 
                "INTS13", 
                "PLAU"),
    downstream = NULL,
    upstream = NULL,
    pAdjustMethod = "bonferroni",
    useModel = "poisson",
    pAdjustThreshold = 0.05,
    logfcThreshold = 0.1)
```

Users can also use `upstream` and `downstream` to specify SNPs proximal to the
gene in the genome.

**Example**

```{r callQTL3_matrix, message=FALSE} 
eqtl3_matrix <- callQTL(
    eQTLObject = eqtl_matrix,
    gene_ids = NULL,
    downstream = -9e7,
    upstream = 2e8,
    pAdjustMethod = "bonferroni",
    useModel = "poisson",
    pAdjustThreshold = 0.05,
    logfcThreshold = 0.05)
```


## Visualize the result.
Here we use `visualizeQTL()` to visualize the result. There are four types of
plots available to visualize sc-eQTL results. Users can choose "histplot",
"violin", "boxplot", or "QTLplot". 

**Example**

```{r visualizeQTL_matrix, message=FALSE}
visualizeQTL(
    eQTLObject = eqtl1_matrix,
    SNPid = "1:632647",
    Geneid = "RPS27",
    groupName = NULL,
    plottype = "QTLplot",
    removeoutlier = TRUE)
```

```{r visualizeQTL_seuratobject, message=FALSE}
visualizeQTL(
    eQTLObject = eqtl1_seurat,
    SNPid = "1:632647",
    Geneid = "RPS27",
    groupName = NULL,
    plottype = "QTLplot",
    removeoutlier = TRUE)
```

In addition, the parameter `groupName` is used to specify a particular
single-cell group of interest.

```{r visualizeQTL_seuratobject_groupName, message=FALSE}
visualizeQTL(
    eQTLObject = eqtl1_seurat,
    SNPid = "1:632647",
    Geneid = "RPS27",
    groupName = "GMP",
    plottype = "QTLplot",
    removeoutlier = TRUE)
```
