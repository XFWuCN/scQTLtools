#' Test Gene Expression Dataset
#'
#' A dataset containing example gene expression data for testing purposes.
#' 100 rows and 2705 colnums. The row names represent gene IDs or SYMBOL and
#' the column names represent cell IDs.
#' @docType data
#' @usage data(testGene)
#' @format A simple \code{matrix}.
#' @keywords datasets
#' @examples
#' data(testGene)
"testGene"

#' Test Genotype Dataset
#'
#' A dataset containing single nucleotide variant data. 1000 rows and 2705
#' colnums. Each row is one variant and each column is one cell.
#' @docType data
#' @usage data(testSNP)
#' @format A simple \code{matrix}.
#' @keywords datasets
#' @examples
#' data(testSNP)
"testSNP"

#' Test SeuratObject
#'
#' A Seurat object for single-cell RNA-seq data.
#' @docType data
#' @usage data(testSeurat)
#' @format A \code{object}
#' @keywords datasets
#' @examples
#' data(testSeurat)
"testSeurat"

#' Test Genotype Dataset
#'
#' A dataset containing single nucleotide variant data.500 rows and 500 colnums
#' Each row is one variant and each column is one cell.
#' @docType data
#' @usage data(testSNP2)
#' @format A simple \code{matrix}.
#' @keywords datasets
#' @examples
#' data(testSNP2)
"testSNP2"

#' Test eqtl object
#'
#' An `eqtlObject` created by the `createQTLObject` function, where the raw
#' expression matrix is normalized using `normalizeGene()`, and both the
#' genotype matrix and the normalized gene expression matrix are filtered
#' by `filterGeneSNP()`.
#' @docType data
#' @usage data(testEQTL)
#' @format A simple \code{object}.
#' @keywords datasets
#' @examples
#' data(testEQTL)
"testEQTL"
