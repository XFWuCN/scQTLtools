#' Test Gene Expression Dataset
#'
#' A dataset containing example gene expression data for testing purposes.
#' It contains 100 rows and 2705 columns. The row names represent gene IDs or
#' SYMBOLs, and the column names represent cell IDs.
#'
#' @name DataSet
#' @aliases testGene
#' @docType data
#' @keywords datasets
#' @format A numeric matrix with 100 rows (genes) and 2705 columns (cells).
NULL

#' Test Genotype Dataset
#'
#' A dataset containing single nucleotide variant data. Each row is one variant
#' and each column is one cell.
#'
#' @name DataSet
#' @aliases testSNP
#' @docType data
#' @keywords datasets
#' @format A numeric matrix with variants as rows and cells as columns.
NULL

#' Test Seurat Object
#'
#' A Seurat object representing single-cell RNA-seq data.
#'
#' @name DataSet
#' @aliases testSeurat
#' @docType data
#' @keywords datasets
#' @format An object of class \code{eQTLObject}.
NULL

#' Test Genotype Dataset
#'
#' A smaller dataset containing single nucleotide variant data for testing.
#' Each row is one variant and each column is one cell.
#'
#' @name DataSet
#' @aliases testSNP2
#' @docType data
#' @keywords datasets
#' @format A numeric matrix with variants as rows and cells as columns.
NULL

#' Test eqtl object
#'
#' An \code{eQTLObject} created using \code{createQTLObject()}, where the raw
#' expression matrix is normalized using \code{normalizeGene()}, and both the
#' genotype matrix and normalized expression matrix are filtered using
#' \code{filterGeneSNP()}.
#'
#' @name DataSet
#' @aliases testEQTL
#' @docType data
#' @keywords datasets
#' An object of class \code{eQTLObject}.
NULL
