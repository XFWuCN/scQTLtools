#' Adjust p-values and filter gene-SNP pairs based on the adjusted p-values.
#'
#' @param result  A data frame that contains information of gene-SNP pairs.
#' @param pAdjustThreshold  Only gene-SNP pairs with adjusted p-values meeting
#' the threshold will be displayed. Default by 0.05.
#' @param pAdjustMethod  Method for p-value adjustment. One of "bonferroni",
#' "holm", "hochberg", "hommel" or "BH". The default option is "bonferroni".
#'
#' @importFrom stats p.adjust
#' @return A data frame with adjusted p-values, filtered by threshold,
#' containing information on gene-SNP pairs.
#' @export
#' @examples
#' example_data <- data.frame(
#'   gene = c("Gene1", "Gene2", "Gene3", "Gene4"),
#'   SNP = c("SNP1", "SNP2", "SNP3", "SNP4"),
#'   pvalue = c(0.001, 0.04, 0.03, 0.0005))
#' pAdjustMethod <- "BH"
#' pAdjustThreshold <- 0.05
#' adjusted_result <- adjust_pvalues(example_data, pAdjustMethod,
#' pAdjustThreshold)
adjust_pvalues <- function(result, pAdjustMethod = "bonferroni",
                            pAdjustThreshold = 0.05) {
    if (!pAdjustMethod %in% c(
        "bonferroni",
        "holm",
        "hochberg",
        "hommel",
        "BH")) {
    stop("Invalid p-adjusted method. Please choose from 'bonferroni', 'holm',
    'hochberg', 'hommel', or'fdr or BH'.")
        }

    result[, "adjusted_pvalue"] <- p.adjust(result[, "pvalue"],
                                            method = pAdjustMethod)
    result <- result[order(result[, "adjusted_pvalue", drop = TRUE]), ,
                    drop = FALSE]
    rownames(result) <- NULL
    result <- result[result$adjusted_pvalue <=
                        pAdjustThreshold, , drop = FALSE]
    result <- result[!is.na(result$adjusted_pvalue), , drop = FALSE]
    return(result)
}
