#' Extract Counts from an Expression Matrix
#'
#' This function retrieves expression counts for a specified gene from an
#' expression matrix, based on the provided list of cells.
#'
#' @param expressionMatrix A matrix containing gene expression data where
#'                         rows represent genes and columns represent cells.
#' @param Geneid A character string or numeric index representing the
#'                specific gene of interest in the expression matrix.
#' @param cells A character vector of cell names (column names of the
#'              expression matrix) from which to extract counts for the
#'              specified gene.
#' @return A numeric vector of expression counts for the specified gene
#'         in the selected cells.
#' @export
#' @examples
#' data(testGene)
#' get_counts(testGene, "CNN2",
#'           c("CGGCAGTGTAGCCCTG", "GGAGGATTCCCGTTCA"))
get_counts <- function(expressionMatrix, Geneid, cells) {
    return(unlist(expressionMatrix[Geneid, cells, drop = TRUE]))
}

#' Retrieve Cells by SNP Value
#'
#' This function extracts the names of cells from a SNP matrix that correspond
#' to a specified value for a given SNP.
#'
#' @param snpMatrix A matrix containing SNP data where rows represent SNPs
#'                  and columns represent cells.
#' @param SNPid A character string or numeric index representing the specific
#'               SNP of interest in the SNP matrix.
#' @param biClassify The user chooses whether to convert the counting method of
#' the snpMatrix to 0/1/2, TRUE indicates conversion,
#' and FALSE indicates no conversion, default is no conversion.
#'
#' @return A list of cell names (column names of the SNP matrix)
#'         that correspond to the specified genotype value for the given SNP.
#' @export
#' @examples
#' data(testSNP)
#' biClassify <- FALSE
#' get_cell_groups(testSNP, "1:632445", biClassify)
get_cell_groups <- function(snpMatrix, SNPid, biClassify) {
    if (biClassify) {
    snpMatrix[snpMatrix == 3] <- 2
    return(list(colnames(snpMatrix)[snpMatrix[SNPid, ] == 1],
            colnames(snpMatrix)[snpMatrix[SNPid, ] == 2]))
    } else {
    return(list(colnames(snpMatrix)[snpMatrix[SNPid, ] == 1],
                colnames(snpMatrix)[snpMatrix[SNPid, ] == 2],
                colnames(snpMatrix)[snpMatrix[SNPid, ] == 3]))
    }
}
