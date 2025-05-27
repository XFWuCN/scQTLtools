#' The package includes five examplt datasets: a gene expression matrix, two SNP
#' genotype matrices, a Seurat object, and an eQTL object.
#'
#' GeneData: A dataset containing example gene expression data. Each row
#' represents a gene and each column represents a cell.
#'
#' @name DataSet
#' @aliases GeneData
#' @docType data
#' @keywords datasets
#' @format A numeric matrix with 100 rows (genes) and 2705 columns (cells).
NULL

#' SNP Genotype Dataset
#'
#' SNPData: A dataset containing single nucleotide (SNP) data. Each row
#' represents a variant and each column represents a cell.
#'
#' @name DataSet
#' @aliases SNPData
#' @docType data
#' @keywords datasets
#' @format A numeric matrix with 1000 rows (SNPs) and 2705 columns (cells).
NULL

#' Seurat Object
#'
#' Seurat_obj: A Seurat object containing \code{GeneData} along with associated
#' metadata stored in the \code{meta.data} slot.
#'
#' @name DataSet
#' @aliases Seurat_obj
#' @docType data
#' @keywords datasets
#' @format An object of class \code{Seurat}.
NULL

#' Genotype Dataset
#'
#' SNPData2: A smaller SNP dataset containing fewer variants and cells. Each row
#' represents a variant and each column represents a cell.
#'
#' @name DataSet
#' @aliases SNPData2
#' @docType data
#' @keywords datasets
#' @format A numeric matrix with 500 rows (SNPs) and 500 columns (cells).
NULL

#' EQTL Object
#'
#' EQTL_obj: An \code{eQTLObject} created using
#' \code{\link{createQTLObject}}, where the raw expression matrix is
#' normalized using \code{\link{normalizeGene}}, and both the genotype matrix
#' and normalized expression matrix were filtered using
#' \code{\link{filterGeneSNP}}.
#'
#' @name DataSet
#' @aliases EQTL_obj
#' @docType data
#' @keywords datasets
#'
#' @seealso
#' \code{\link{createQTLObject}}, \code{\link{normalizeGene}} and
#' \code{\link{filterGeneSNP}}, for the underlying functions that do the work.
#'
#' An object of class \code{eQTLObject}.
NULL
