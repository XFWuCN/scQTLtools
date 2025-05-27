#' Normalize gene expression data.
#'
#' This function performs normalization of the gene expression matrix stored
#' within an \code{eQTLObject}, aiming to remove technical biases and
#' improve comparability across cells. Multiple normalization methods are
#' supported, including "logNormalize", "CPM", "TPM", "DESeq", and "limma".
#'
#' @param eQTLObject  An S4 object of class \code{eQTLObject} containing the raw
#' gene expression matrix.
#' @param method  Character string specifying the normalization method to use.
#' Must be one of \code{"logNormalize"}, \code{"CPM"}, \code{"TPM"},
#' \code{"DESeq"}, or \code{"limma"}. Default is \code{"logNormalize"}.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @importFrom limma normalizeBetweenArrays
#'
#' @return An \code{eQTLObject} with the normalized gene expression matrix
#' stored in the slot named \code{"normExpMat"}.
#'
#' @export
#' @examples
#' data(EQTL_obj)
#' eqtl <- normalizeGene(EQTL_obj, method = "logNormalize")
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
