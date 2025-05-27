#' Fit Linear Model for eQTL Mapping
#'
#' This function performs linear regression to identify gene–SNP associations
#' based on single-cell expression and genotype data stored in an
#' \code{eQTLObject}.
#' @param eQTLObject An S4 object of class \code{eQTLObject}.
#' @param geneIDs Character vector of gene IDs to include in the model fitting.
#' @param snpIDs Character vector of SNP IDs to include in the model fitting.
#' @param biClassify Logical; whether to convert genotype encoding in snpMatrix
#' to 0, 1, and 2. \code{TRUE} indicates conversion; \code{FALSE} indicates no
#' conversion (default).
#' @param pAdjustThreshold  Only SNP–gene pairs with adjusted p-values below
#' the threshold will be retained. Default is 0.05.
#' @param pAdjustMethod  Method used for multiple testing correction. One of
#' \code{"bonferroni"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, or
#' \code{"BH"}. Default is \code{"bonferroni"}.
#' @param logfcThreshold  The minimum beta coefficient (effect size) required
#' to report a SNP–gene pair as an eQTL.
#'
#' @importFrom stats lm
#' @return A data frame of gene–SNP pairs that pass the filtering criteria,
#' including beta coefficients, p-values, adjusted p-values, and group labels.
#' @export
#' @examples
#' data(EQTL_obj)
#' Gene <- rownames(slot(EQTL_obj, "filterData")$expMat)
#' SNP <- rownames(slot(EQTL_obj, "filterData")$snpMat)
#' linearResult <- linearModel(
#'   eQTLObject = EQTL_obj,
#'   geneIDs = Gene,
#'   snpIDs = SNP,
#'   biClassify = FALSE,
#'   pAdjustMethod = "bonferroni",
#'   pAdjustThreshold = 0.05,
#'   logfcThreshold = 0.025)
linearModel <- function(eQTLObject, geneIDs, snpIDs,
    biClassify = FALSE, pAdjustMethod = "bonferroni",
    pAdjustThreshold = 0.05, logfcThreshold = 0.1) {
    unique_group <- unique(load_group_info(eQTLObject)[["group"]])

    message("Start the sc-eQTLs calling.")

    result_all <- do.call(rbind, lapply(unique_group, function(k) {
        process_group_common(eQTLObject, geneIDs, snpIDs, biClassify, k,
                            model = "linear")
    }))

    message("Finished!")

    result_all <- adjust_pvalues(result_all, pAdjustMethod, pAdjustThreshold)
    result_all <- filter_by_abs_b(result_all, logfcThreshold)

    return(result_all)
}
