#' Poisson model fitting the gene expression matrix and genotype matrix.
#'
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param geneIDs Matching genes can be used to fit data.
#' @param snpIDs Matching SNPs can be used to fit data.
#' @param biClassify The user chooses whether to convert the counting method of
#' the snpMatrix to 0/1/2, TRUE indicates conversion, and FALSE indicates no
#' conversion, default is FALSE.
#' @param pAdjustMethod Methods for p-value adjusting, one of "bonferroni",
#' "holm", "hochberg", "hommel" or "BH". The default option is "bonferroni".
#' @param pAdjustThreshold  Only SNP-Gene pairs with adjusted p-values meeting
#' the threshold will be displayed. The default value is 0.05.
#' @param logfcThreshold Represents the minimum beta threshold for fitting
#' SNP-Gene pairs.
#'
#' @importFrom stats p.adjust poisson glm
#' @return Dataframe that contains gene-SNP pairs' information.
#' @export
#' @examples
#' data(testEQTL)
#' Gene <- rownames(slot(testEQTL, "filterData")$expMat)
#' SNP <- rownames(slot(testEQTL, "filterData")$snpMat)
#' poissonResult <- poissonModel(
#'   eQTLObject = testEQTL,
#'   geneIDs = Gene,
#'   snpIDs = SNP,
#'   biClassify = FALSE,
#'   pAdjustMethod = "bonferroni",
#'   pAdjustThreshold = 0.05,
#'   logfcThreshold = 0.025
#' )
poissonModel <- function(eQTLObject, geneIDs, snpIDs,
    biClassify = FALSE, pAdjustMethod = "bonferroni",
    pAdjustThreshold = 0.05, logfcThreshold = 0.1) {
    unique_group <- unique(load_group_info(eQTLObject)[["group"]])

    message("Start the sc-eQTLs calling.")

    result_all <- do.call(rbind, lapply(unique_group,
        function(k) {
            process_group_common(eQTLObject, geneIDs, snpIDs, biClassify, k,
                                model = "poisson")
        }))

    message("Finished!")

    result_all <- adjust_pvalues(result_all, pAdjustMethod, pAdjustThreshold)
    result_all <- filter_by_abs_b(result_all, logfcThreshold)

    return(result_all)
}

