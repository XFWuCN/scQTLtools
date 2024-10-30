#' normalizeGene: Normalize the gene expression data.
#' @description Gene expression matrix normalization is necessary to eliminate
#' technical biases and variabilities, ensuring accurate and comparable
#' analysis of gene expression data. Here we provide `normalizeGene()`to
#' normalize the data.
#'
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param method Method for normalizing for gene expression dataframe, one of
#' "logNormalize", "CPM", "TPM", "DESeq" or "limma"
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @importFrom limma normalizeBetweenArrays
#'
#' @return A normalized gene expression matrix.
#'
#' @export
#' @examples
#' data(testEQTL)
#' eqtl <- normalizeGene(testEQTL, method = "logNormalize")
normalizeGene <- function(eQTLObject, method = "logNormalize") {
    expressionMatrix <- get_raw_data(eQTLObject)[["rawExpMat"]]
    method <- method
    if (!method %in% c("logNormalize", "CPM", "TPM", "DESeq", "limma")) {
        stop("Invalid method.
        Please choose from 'logNormalize', 'CPM, 'TPM', 'DESeq' or 'limma' .")
    }
    rowsum <- apply(expressionMatrix, 1, sum)
    expressionMatrix <- expressionMatrix[rowsum != 0, , drop = FALSE]

    normalizedData <- switch(method,
                            logNormalize = log_normalize(expressionMatrix),
                            CPM = CPM_normalize(expressionMatrix),
                            TPM = TPM_normalize(expressionMatrix),
                            DESeq = DESeq_normalize(expressionMatrix),
                            limma = limma_normalize(expressionMatrix))

    message("Normalization completed using method: ", method, "\n")
    message("Dimensions of normalized data:",
            paste(dim(normalizedData),
                    collapse = " "),
            "\n")
    eQTLObject <- set_raw_data(eQTLObject, normalizedData, "normExpMat")
    return(eQTLObject)
}
