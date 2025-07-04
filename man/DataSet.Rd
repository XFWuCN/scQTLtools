% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{DataSet}
\alias{DataSet}
\alias{GeneData}
\alias{SNPData}
\alias{Seurat_obj}
\alias{SNPData2}
\alias{EQTL_obj}
\title{The package includes five examplt datasets: a gene expression matrix, two SNP
genotype matrices, a Seurat object, and an eQTL object.}
\format{
A numeric matrix with 100 rows (genes) and 2705 columns (cells).

A numeric matrix with 1000 rows (SNPs) and 2705 columns (cells).

An object of class \code{Seurat}.

A numeric matrix with 500 rows (SNPs) and 500 columns (cells).
}
\description{
GeneData: A dataset containing example gene expression data. Each row
represents a gene and each column represents a cell.

SNPData: A dataset containing single nucleotide (SNP) data. Each row
represents a variant and each column represents a cell.

Seurat_obj: A Seurat object containing \code{GeneData} along with associated
metadata stored in the \code{meta.data} slot.

SNPData2: A smaller SNP dataset containing fewer variants and cells. Each row
represents a variant and each column represents a cell.

EQTL_obj: An \code{eQTLObject} created using
\code{\link{createQTLObject}}, where the raw expression matrix is
normalized using \code{\link{normalizeGene}}, and both the genotype matrix
and normalized expression matrix were filtered using
\code{\link{filterGeneSNP}}.
}
\seealso{
\code{\link{createQTLObject}}, \code{\link{normalizeGene}} and
\code{\link{filterGeneSNP}}, for the underlying functions that do the work.

An object of class \code{eQTLObject}.
}
\keyword{datasets}
