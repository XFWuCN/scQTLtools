#' Normalize the gene expression matrix with DESeq.
#'
#' `DESeq_normalize()`normalizes an expression matrix using the DESeq2 package.
#'
#' @param expressionMatrix Input raw gene expression matrix.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @return A gene expression matrix after normalized.
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
    dds <- DESeq2::DESeq(dds)
    normalizedData <- DESeq2::counts(dds, normalized = TRUE)
    return(normalizedData)
}
