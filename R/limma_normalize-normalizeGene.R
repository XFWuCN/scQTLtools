#' Normalize the gene expression matrix with limma
#'
#' \code{limma_normalize()} normalizes an expression matrix using the quantile
#' normalization method provided by the limma package.
#'
#' @param expressionMatrix  A numeric matrix of raw gene expression counts,
#' with genes as rows and cells as columns.
#' @importFrom limma normalizeBetweenArrays
#' @return A normalized gene expression matrix after applying limma
#' normalization.
#' @export
#'
#' @examples
#' data(GeneData)
#' limma_normalize(GeneData)
limma_normalize <- function(expressionMatrix) {
    normalizedData <- normalizeBetweenArrays(expressionMatrix,
                                                    method = "quantile")
    return(normalizedData)
}
