#' Normalize gene expression using CPM.
#'
#' \code{CPM_normalize()} applies Counts Per Million (CPM) normalization to a
#' raw gene expression matrix.
#'
#' @param expressionMatrix  A numeric matrix of raw gene expression counts,
#' with genes as rows and cells as columns.
#'
#' @return A normalized gene expression matrix after applying CPM normalization.
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
