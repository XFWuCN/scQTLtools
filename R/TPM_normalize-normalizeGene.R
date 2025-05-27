#' Normalize the gene expression matrix with TPM
#'
#' \code{TPM_normalize()} scales an expression matrix using Transcripts Per
#' Million (TPM) normalization, applying logarithm and scaling operations to
#' adjust data based on library size.
#'
#' @param expressionMatrix  A numeric matrix of raw gene expression counts,
#' with genes as rows and cells as columns.
#'
#' @return A normalized gene expression matrix after applying TPM
#' normalization.
#' @export
#'
#' @examples
#' data(GeneData)
#' TPM_normalize(GeneData)
TPM_normalize <- function(expressionMatrix) {
    library_size <- colSums(expressionMatrix)
    normalizedData <- log1p(expressionMatrix / library_size * 1e6)
    return(normalizedData)
}
