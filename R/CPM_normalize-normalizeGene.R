#' Normalize the gene expression matrix with CPM.
#'
#'`CPM_normalize()` scales an expression matrix using Counts Per Million (CPM)
#' normalization, applying logarithm and scaling operations to adjust data.
#'
#' @param expressionMatrix Input raw gene expression matrix.
#'
#' @return A gene expression matrix after normalized.
#' @export
#'
#' @examples
#' data(testGene)
#' CPM_normalize(testGene)
CPM_normalize <- function(expressionMatrix) {
    total_counts <- colSums(expressionMatrix)
    normalizedData <- log1p(
        sweep(expressionMatrix, 2, total_counts, "/") * 1e6)
    return(normalizedData)
}
