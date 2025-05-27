#' Normalize the gene expression matrix with logNormalize method.
#'
#'\code{log_normalize()} transforms an expression matrix by applying logarithm
#' and scaling operations to normalize data.
#'
#' @param expressionMatrix  A numeric matrix of raw gene expression counts,
#' with genes as rows and cells as columns.
#'
#' @return A normalized gene expression matrix after applying logNormalize
#' normalization.
#' @export
#'
#' @examples
#' data(GeneData)
#' log_normalize(GeneData)
log_normalize <- function(expressionMatrix) {
    normalizedData <- log1p(sweep(expressionMatrix,
                                2,
                                colSums(expressionMatrix),
                                FUN = "/") * 1e4)
    return(normalizedData)
}
