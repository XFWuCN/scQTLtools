#' Normalize the gene expression matrix with DESeq2.
#'
#' \code{DESeq_normalize()} normalizes a raw gene expression matrix using the
#' \pkg{DESeq2} package.
#'
#' @param expressionMatrix  A numeric matrix of raw gene expression counts,
#' with genes as rows and cells as columns.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @return A normalized gene expression matrix after applying DESeq
#' normalization.
#' @export
#'
#' @examples
#' data(testGene)
#' DESeq_normalize(testGene)
DESeq_normalize <- function(expressionMatrix) {
    sampleDataframe <- colnames(expressionMatrix)
    sampleDataframe <- as.data.frame(sampleDataframe)
    options(warn = -1)
    dds <- DESeqDataSetFromMatrix(countData = expressionMatrix,
                                colData = sampleDataframe,
                                design = ~1)
    dds <- DESeq(dds)
    normalizedData <- counts(dds, normalized = TRUE)
    return(normalizedData)
}
