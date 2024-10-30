#' filterGeneSNP: Filter gene expression matrix and genotype matrix.
#'
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param snpNumOfCellsPercent Only SNPs where cells with each of the
#' different genotypes (REF and ALT, or AA, Aa, and aa) individually account
#' for at least `snpNumOfCellsPercent`% of the total cells will be considered.
#' Default by 10.
#' @param expressionMin threshold for valid gene expression levels, utilized
#' alongside another parameter,expression.number.of.cells.Default by 0.
#' @param expressionNumOfCellsPercent Only genes with expression levels
#' exceeding `expressionMin` in at least `expressionNumOfCellsPercent`%
#' of cells are considered. The default value is 10.
#' @importFrom progress progress_bar
#'
#' @return filtered matrices.
#' @export
#'
#' @examples
#' data(testSNP)
#' data(testGene)
#' eqtl <- createQTLObject(snpMatrix = testSNP, genedata = testGene)
#' eqtl <- normalizeGene(eqtl)
#' eqtl <- filterGeneSNP(eqtl,
#'   snpNumOfCellsPercent = 2,
#'   expressionMin = 0,
#'   expressionNumOfCellsPercent = 2)
filterGeneSNP <- function(eQTLObject,
                        snpNumOfCellsPercent = 2,
                        expressionMin = 0,
                        expressionNumOfCellsPercent = 2) {
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
