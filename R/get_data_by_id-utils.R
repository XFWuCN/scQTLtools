#' Extract Counts from an Expression Matrix
#'
#' This function retrieves expression counts for a specified gene from an
#' expression matrix, based on the provided list of cells.
#'
#' @param expressionMatrix A numeric matrix of gene expression counts,
#' with genes as rows and cells as columns.
#' @param Geneid A character string or numeric index representing the specific
#' gene of interest in \code{expressionMatrix}.
#' @param cells A character vector of cell names (column names of
#' \code{expressionMatrix}) from which to extract counts for the specified gene.
#'
#' @return A numeric vector containing the expression counts of the specified
#' gene in the selected cells.
#' @export
#' @examples
#' data(GeneData)
#' get_counts(GeneData, "CNN2",
#'           c("CGGCAGTGTAGCCCTG", "GGAGGATTCCCGTTCA"))
get_counts <- function(expressionMatrix, Geneid, cells) {
    return(unlist(expressionMatrix[Geneid, cells, drop = TRUE]))
}

#' Retrieve Cells by SNP Value
#'
#' This function extracts the names of cells from a SNP matrix that correspond
#' to a specified value for a given SNP.
#'
#' @param snpMatrix A genotype matrix where each row is a snp and each column
#' is a cell. Encoding should be 0, 1, 2, 3.
#' @param SNPid A character string or numeric index representing the specific
#'               SNP of interest in the SNP matrix.
#' @param biClassify Logical; whether to convert genotype encoding in snpMatrix
#' to 0, 1, and 2. \code{TRUE} indicates conversion; \code{FALSE} indicates no
#' conversion (default).
#'
#' @return A list of character vectors. Each vector contains the names of cells
#' (i.e., column names of \code{snpMatrix}) corresponding to a specific
#' genotype value at the given SNP.
#' @export
#' @examples
#' data(SNPData)
#' biClassify <- FALSE
#' get_cell_groups(SNPData, "1:632445", biClassify)
get_cell_groups <- function(snpMatrix, SNPid, biClassify) {
    if (biClassify) {
    snpMatrix[snpMatrix == 3] <- 2
    return(list(colnames(snpMatrix)[snpMatrix[SNPid, ] == 1],
            colnames(snpMatrix)[snpMatrix[SNPid, ] == 2]))
    } else {
    return(list(colnames(snpMatrix)[snpMatrix[SNPid, ] == 1],
                colnames(snpMatrix)[snpMatrix[SNPid, ] == 3],
                colnames(snpMatrix)[snpMatrix[SNPid, ] == 2]))
    }
}
