% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_normalize-normalizeGene.R
\name{log_normalize}
\alias{log_normalize}
\title{Normalize the gene expression matrix with logNormalize method.}
\usage{
log_normalize(expressionMatrix)
}
\arguments{
\item{expressionMatrix}{A numeric matrix of raw gene expression counts,
with genes as rows and cells as columns.}
}
\value{
A normalized gene expression matrix after applying logNormalize
normalization.
}
\description{
\code{log_normalize()} transforms an expression matrix by applying logarithm
and scaling operations to normalize data.
}
\examples{
data(GeneData)
log_normalize(GeneData)
}
