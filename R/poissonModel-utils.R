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
poissonModel <- function(
    eQTLObject,
    geneIDs,
    snpIDs,
    biClassify = FALSE,
    pAdjustMethod = "bonferroni",
    pAdjustThreshold = 0.05,
    logfcThreshold = 0.1) {

    expressionMatrix <- round(get_filter_data(eQTLObject)[["expMat"]] * 1000)
    snpMatrix <- get_filter_data(eQTLObject)[["snpMat"]]
    unique_group <- unique(load_group_info(eQTLObject)[["group"]])

    result_all <- data.frame()

    message("Start the sc-eQTLs calling.")
    result_all <- do.call(rbind, lapply(unique_group, function(k) {
    result <- data.frame(
        SNPid = character(),
        group = character(),
        Geneid = character(),
        pvalue = double(),
        adjusted_pvalue = double(),
        b = double(),
        abs_b = double(),
        Remark = character(),
        stringsAsFactors = FALSE)

    split_cells <- rownames(load_group_info(eQTLObject)
                            )[load_group_info(eQTLObject)[["group"]] == k]
    expressionMatrix_split <- expressionMatrix[, split_cells, drop = FALSE]
    snpMatrix_split <- snpMatrix[, split_cells, drop = FALSE]

    replace_2_and_3 <- function(x) {
        ifelse(x == 2, 3, ifelse(x == 3, 2, x))
    }

    if (biClassify == FALSE) {
        snpMatrix_split <- as.data.frame(snpMatrix_split)
        snpMatrix_split <- snpMatrix_split %>%
            mutate_all(list(~ replace_2_and_3(.)))
        snpMatrix_split <- as.matrix(snpMatrix_split)
    } else {
    snpMatrix_split[snpMatrix_split == 3] <- 2
    }

    pb <- initialize_progress_bar(length(snpIDs), k)

    result_list <- lapply(seq_len(length(snpIDs)), function(i) {
        snp_mat <- process_matrix(snpIDs[i], snpMatrix_split, "snp_mat")

        gene_results <- lapply(seq_len(length(geneIDs)), function(j) {
            gene_mat <- process_matrix(geneIDs[j], expressionMatrix_split,
                                        "gene_mat")

            combined_df <- merge(snp_mat, gene_mat, by = "cells")
            combined_df <- subset(combined_df, snp_mat != 0)

            lmodel <- glm(combined_df$gene_mat ~ combined_df$snp_mat,
                        family = poisson())

            if(length(summary(lmodel)$coefficients[, "Pr(>|z|)"]) >= 2) {
                lmout_pvalue <- summary(lmodel)$coefficients[2, "Pr(>|z|)"]
                lmout_b <- summary(lmodel)$coefficients[2, "Estimate"]
                new_row <- data.frame(
                    SNPid = snpIDs[i],
                    group = k,
                    Geneid = geneIDs[j],
                    pvalue = lmout_pvalue,
                    b = lmout_b)
            return(new_row)
            } else {
                return(NULL)
            }
        })
        pb$tick()
        gene_results <- do.call(rbind, gene_results)
        gene_results <- gene_results[vapply(gene_results,
                                            is.null,
                                            logical(1)) == FALSE, ]
        return(gene_results)
})
    result <- do.call(rbind, result_list)
    message("Finished!")

    result <- adjust_pvalues(result, pAdjustThreshold, pAdjustMethod)
    result <- filter_by_abs_b(result, logfcThreshold)
    }))
    return(result_all)
}

