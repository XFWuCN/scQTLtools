#' Normalize the gene expression matrix with limma
#'
#' `limma_normalize()`normalizes an expression matrix using the quantile
#' normalization method provided by the limma package.
#'
#' @param expressionMatrix Input raw gene expression matrix.
#' @importFrom limma normalizeBetweenArrays
#' @return A gene expression matrix after normalized.
#' @export
#'
#' @examples
#' data(testGene)
#' limma_normalize(testGene)
limma_normalize <- function(expressionMatrix) {
    normalizedData <- normalizeBetweenArrays(expressionMatrix,
                                                    method = "quantile")
    return(normalizedData)
}
