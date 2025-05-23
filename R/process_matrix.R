#' Process a matrix to extract a row and convert it to a data frame.
#'
#' @param id Character string specifying the row name to extract from the
#' matrix.
#' @param matrix A matrix from which the row will be extracted.
#' @param name Character string specifying the column name for the extracted
#' values.
#'
#' @return A data frame containing the extracted row and a column with the
#' row names.
#' @export
#'
#' @examples
#' rownames <- c("CNN2", "TIGD2", "DTD2")
#' colnames <- c("Col1", "Col2", "Col3", "Col4")
#' matrix_data <- matrix(1:12, nrow = 3, ncol = 4,
#'   dimnames = list(rownames, colnames))
#' geneid <- "CNN2"
#' gene_mat <- process_matrix(geneid, matrix_data, "gene_mat")
process_matrix <- function(id, matrix, name){
    data <- matrix[id, , drop = TRUE]
    data <- as.data.frame(data)
    colnames(data) <- name
    data$cells <- rownames(data)
    return(data)
}
