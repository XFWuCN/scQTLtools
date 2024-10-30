#' Normalize the gene expression matrix with logNormalize method.
#'
#'`log_normalize()` transforms an expression matrix by applying logarithm and
#' scaling operations to normalize data.
#'
#' @param expressionMatrix Input raw gene expression matrix.
#'
#' @return A gene expression matrix after normalized.
#' @export
#'
#' @examples
#' data(testGene)
#' log_normalize(testGene)
log_normalize <- function(expressionMatrix) {
    normalizedData <- log1p(sweep(expressionMatrix,
                                2,
                                colSums(expressionMatrix),
                                FUN = "/") * 1e4)
    return(normalizedData)
}
