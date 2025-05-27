#' Filter gene expression and genotype matrices by cell percentage thresholds.
#'
#' @param eQTLObject An S4 object of class \code{eQTLObject}.
#' @param snpNumOfCellsPercent Numeric. Minimum percentage of cells required
#' for each SNP genotype (e.g., AA, AG, and GG). Only SNPs where each genotype
#' occurs in at least this proportion of cells are retained. Default is 10.
#' @param expressionMin Numeric. Expression threshold used in combination with
#' \code{expressionNumOfCellsPercent} to filter lowly expressed genes.
#' Default is 0.
#' @param expressionNumOfCellsPercent Numeric. Minimum percentage of cells in
#' which a gene's expression must exceed \code{expressionMin} for the gene to
#' be retained. Default is 10.
#' @importFrom progress progress_bar
#'
#' @return An updated \code{eQTLObject} with filtered gene expression and SNP
#' matrices.
#' @export
#'
#' @examples
#' data(SNPData)
#' data(GeneData)
#' eqtl <- createQTLObject(snpMatrix = SNPData, genedata = GeneData)
#' eqtl <- normalizeGene(eqtl)
#' eqtl <- filterGeneSNP(eqtl,
#'   snpNumOfCellsPercent = 2,
#'   expressionMin = 0,
#'   expressionNumOfCellsPercent = 2)
filterGeneSNP <- function(eQTLObject,
                        snpNumOfCellsPercent = 10,
                        expressionMin = 0,
                        expressionNumOfCellsPercent = 10) {
    filtered_expressionMatrix <- filter_expr(eQTLObject = eQTLObject,
                    expressionMin = expressionMin,
                    expressionNumOfCellsPercent = expressionNumOfCellsPercent)

    snpMatrix <- as.data.frame(get_raw_data(eQTLObject)[["snpMat"]])
    snp.list <- rownames(snpMatrix)
    snp.number.of.cells <- ceiling(
        (snpNumOfCellsPercent / 100) * ncol(snpMatrix)
        )
    biClassify <- load_biclassify_info(eQTLObject)
    if (biClassify == TRUE) {
        snpMatrix[snpMatrix == 3] <- 2

        snp_counts <- apply(snpMatrix, 1, function(row) {
            c(sum(row == 1), sum(row == 2))
        })
    snp_counts_df <- as.data.frame(t(snp_counts))
    names(snp_counts_df) <- c("count_ref", "count_alt")

    filtered_snp_ids <- rownames(snp_counts_df)[
        apply(snp_counts_df, 1, function(row) {
        all(row > snp.number.of.cells)
        })
    ]
    filtered_snpMatrix <- snpMatrix[filtered_snp_ids, , drop = FALSE]
    } else if (biClassify == FALSE) {
        snp_counts <- apply(snpMatrix, 1, function(row) {
            c(sum(row == 1), sum(row == 3), sum(row == 2))
    })
    snp_counts_df <- as.data.frame(t(snp_counts))
    names(snp_counts_df) <- c("count_AA", "count_Aa", "count_aa")

    filtered_snp_ids <- rownames(snp_counts_df)[
        apply(snp_counts_df, 1, function(row) {
            all(row > snp.number.of.cells)
        })
    ]
    filtered_snpMatrix <- snpMatrix[filtered_snp_ids, , drop = FALSE]
    filtered_snpMatrix <- as.matrix(filtered_snpMatrix)
    }

    eQTLObject <- set_filter_data(eQTLObject, filtered_expressionMatrix,
                                    "expMat")
    eQTLObject <- set_filter_data(eQTLObject, filtered_snpMatrix, "snpMat")
    return(eQTLObject)
}


# @rdname filterGeneSNP_internals
filter_expr <- function(eQTLObject,
                        expressionMin,
                        expressionNumOfCellsPercent){
    if (!is.null(get_raw_data(eQTLObject)[["normExpMat"]])) {
        nor_expressionMatrix <-
            as.data.frame(get_raw_data(eQTLObject)[["normExpMat"]])
        expression.number.of.cells <- ceiling(
            (expressionNumOfCellsPercent / 100) * ncol(nor_expressionMatrix)
        )
        valid_genes <- rownames(nor_expressionMatrix)[
            apply(nor_expressionMatrix > expressionMin, 1, sum) >=
            expression.number.of.cells
        ]
        filtered_expressionMatrix <- nor_expressionMatrix[valid_genes, ,
                                                        drop = FALSE]
        filtered_expressionMatrix <- as.matrix(filtered_expressionMatrix)
    } else {
        stop("Please normalize the raw expression data first.")
    }
    return(filtered_expressionMatrix)
}
