#' Normalize the gene expression matrix with TPM
#'
#' `TPM_normalize()` scales an expression matrix using Transcripts Per Million
#' (TPM) normalization, applying logarithm and scaling operations to adjust
#' data based on library size.
#'
#' @param expressionMatrix Input raw gene expression matrix.
#'
#' @return A gene expression matrix after normalized.
#' @export
#'
#' @examples
#' data(testGene)
#' TPM_normalize(testGene)
TPM_normalize <- function(expressionMatrix) {
    library_size <- colSums(expressionMatrix)
    normalizedData <- log1p(expressionMatrix / library_size * 1e6)
    return(normalizedData)
}
