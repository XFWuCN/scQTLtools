#' Remove outliers from gene expression data and update cell groups.
#'
#' This function identifies and removes outlier expression values for a
#' specified gene based on the Median Absolute Deviation (MAD) method.
#' It then filters provided genotype-based cell groups, returning only
#' cells with non-outlier expression values.
#'
#' @param exprsMat Input gene expression matrix with genes as rows and cells as
#' columns.
#' @param Geneid Character string specifying the gene ID to examine.
#' @param A_cells Character vector of cell names belonging to genotype group A.
#' @param B_cells Character vector of cell names belonging to genotype group B.
#' @param C_cells Optional character vector of cell names belonging to genotype
#' group C; if \code{NULL}, function returns two genotype groups.
#' @importFrom stats median
#' @return A named list of filtered cell vectors.
#' @export
#'
#' @examples
#' ## Mock expression matrix
#' set.seed(123)
#' exprsMat <- matrix(rnorm(200), nrow = 5)
#' rownames(exprsMat) <- paste0("Gene", 1:nrow(exprsMat))
#' colnames(exprsMat)  <- paste0("cell", 1:ncol(exprsMat))
#' A_cells <- colnames(exprsMat)[1:13] # Example A cell list
#' B_cells <- colnames(exprsMat)[14:26]  # Example B cell list
#' C_cells <- colnames(exprsMat)[27:40]  # Example C cell list
#' remove_outliers(exprsMat, "Gene1", A_cells, B_cells, C_cells)
remove_outliers <- function(exprsMat,
                            Geneid,
                            A_cells,
                            B_cells,
                            C_cells = NULL) {
    sample_gene <- exprsMat[Geneid, , drop = FALSE]
    sample_no_zero <- sample_gene[sample_gene != 0]
    med <- median(unlist(sample_no_zero))
    mad <- mad(unlist(sample_no_zero))
    filter <- as.data.frame(sample_no_zero[sample_no_zero < med + 4 * mad])
    if (!is.null(C_cells)) {
        AA_cells <- intersect(A_cells, rownames(filter))
        Aa_cells <- intersect(B_cells, rownames(filter))
        aa_cells <- intersect(C_cells, rownames(filter))
        return(list(AA_cells = AA_cells,
                    Aa_cells = Aa_cells,
                    aa_cells = aa_cells))
        } else {
        ref_cells <- intersect(A_cells, rownames(filter))
        alt_cells <- intersect(B_cells, rownames(filter))
        return(list(ref_cells = ref_cells,
                alt_cells = alt_cells))
    }
}
