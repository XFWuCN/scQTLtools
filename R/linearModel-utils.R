#' Linear model fitting the gene expression matrix and genotype matrix.
#'
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param geneIDs Matching genes can be used to fit data.
#' @param snpIDs Matching SNPs can be used to fit data.
#' @param biClassify The user chooses whether to convert the counting method of
#' the snpMatrix to 0/1/2, TRUE indicates conversion, and FALSE indicates no
#' conversion, default is no conversion.
#' @param pAdjustMethod  Methods for p-value adjusting, one of "bonferroni",
#' "holm", "hochberg", "hommel" or "BH". The default option is "bonferroni".
#' @param pAdjustThreshold  Only SNP-Gene pairs with adjusted p-values meeting
#' the threshold will be displayed. Default by 0.05.
#' @param logfcThreshold Represents the minimum beta threshold for fitting
#' SNP-Gene pairs. Default by 0.1.
#'
#' @importFrom stats lm p.adjust
#' @return Dataframe that contains gene-SNP pairs' information.
#' @export
#' @examples
#' data(testEQTL)
#' Gene <- rownames(slot(testEQTL, "filterData")$normExpMat)
#' SNP <- rownames(slot(testEQTL, "filterData")$snpMat)
#' linearResult <- linearModel(
#'   eQTLObject = testEQTL,
#'   geneIDs = Gene,
#'   snpIDs = SNP,
#'   biClassify = FALSE,
#'   pAdjustMethod = "bonferroni",
#'   pAdjustThreshold = 0.05,
#'   logfcThreshold = 0.025)
linearModel <- function(
    eQTLObject,
    geneIDs,
    snpIDs,
    biClassify = FALSE,
    pAdjustMethod = "bonferroni",
    pAdjustThreshold = 0.05,
    logfcThreshold = 0.1) {
    if (length(eQTLObject@filterData) == 0) {
        stop("Please filter the data first.")
    }
    expressionMatrix <- round(eQTLObject@filterData$expMat * 1000)
    snpMatrix <- eQTLObject@filterData$snpMat
    unique_group <- unique(eQTLObject@groupBy$group)

    result_all <- data.frame()

    message("Start calling eQTL")
    for (k in unique_group) {
    result <- data.frame(
        SNPid = character(),
        group = character(),
        Geneid = character(),
        pvalue = double(),
        adjusted_pvalue = double(),
        b = double(),
        abs_b = double(),
        Remark = character(),
        stringsAsFactors = FALSE
    )

    split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == k]
    expressionMatrix_split <- expressionMatrix[, split_cells]
    snpMatrix_split <- snpMatrix[, split_cells]

    if (biClassify == FALSE) {
        message(k, ":")
        message("0%   10   20   30   40   50   60   70   80   90   100%")
        message("[----|----|----|----|----|----|----|----|----|----|")
        pb <- progress_bar$new(
            total = length(snpIDs),
            format = "[:bar]",
            clear = FALSE,
            width = 51)
    for (i in seq_len(length(snpIDs))) {
        snpid <- snpIDs[i]
        snp_mat <- snpMatrix_split[snpid, ]
        snp_mat <- as.data.frame(snp_mat)
        snp_mat$cells <- rownames(snp_mat)

        replace_2_and_3 <- function(x) {
            ifelse(x == 2, 3, ifelse(x == 3, 2, x))
        }

        snp_mat_new <- snp_mat %>%
            mutate_all(list(~ replace_2_and_3(.)))

        for (j in seq_len(length(geneIDs))) {
            gene_id <- geneIDs[j]
            gene_mat <- expressionMatrix_split[gene_id, ]
            gene_mat <- as.data.frame(gene_mat)
            gene_mat$cells <- rownames(gene_mat)

            combined_df <- merge(snp_mat, gene_mat, by = "cells")
            combined_df <- subset(combined_df, snp_mat != 5)

            lmodel <- lm(gene_mat ~ snp_mat, data = combined_df)

            if(length(summary(lmodel)$coefficients[, "Pr(>|t|)"]) >= 2) {
                lmout_pvalue <- summary(lmodel)$coefficients[2, "Pr(>|t|)"]
                lmout_b <- summary(lmodel)$coefficients[2, "Estimate"]
                new_row <- data.frame(
                SNPid = snpid,
                group = k,
                Geneid = geneIDs[j],
                pvalue = lmout_pvalue,
                b = lmout_b)
            result <- rbind(result, new_row)
            }
        }
        pb$tick()
    }
        message("finished!")
    } else if (biClassify == TRUE) {
        snpMatrix_split[snpMatrix_split == 3] <- 2
        message(k, ":")
        message("0%   10   20   30   40   50   60   70   80   90   100%")
        message("[----|----|----|----|----|----|----|----|----|----|")
        pb <- progress_bar$new(
            total = length(snpIDs),
            format = "[:bar]",
            clear = FALSE,
            width = 51)
    for (i in seq_len(length(snpIDs))) {
        snpid <- snpIDs[i]
        snp_mat <- snpMatrix_split[snpid, ]
        snp_mat <- as.data.frame(snp_mat)
        snp_mat$cells <- rownames(snp_mat)

        for (j in seq_len(length(geneIDs))) {
            gene_id <- geneIDs[j]
            gene_mat <- expressionMatrix_split[gene_id, ]
            gene_mat <- as.data.frame(gene_mat)
            gene_mat$cells <- rownames(gene_mat)

            combined_df <- merge(snp_mat, gene_mat, by = "cells")
            combined_df <- subset(combined_df, snp_mat != 5)

            lmodel <- lm(gene_mat ~ snp_mat, data = combined_df)

            if(length(summary(lmodel)$coefficients[, "Pr(>|t|)"]) >= 2) {
            lmout_pvalue <- summary(lmodel)$coefficients[2, "Pr(>|t|)"]
            lmout_b <- summary(lmodel)$coefficients[2, "Estimate"]
            new_row <- data.frame(
                SNPid = snpid,
                group = k,
                Geneid = geneIDs[j],
                pvalue = lmout_pvalue,
                b = lmout_b)
            result <- rbind(result, new_row)
            }
        }
        pb$tick()
        }
        message("finished!")
    } else {
    stop("biClassify can only be selected as 'TRUE' or 'FALSE'")
    }

    # p-value correction methods
    if (!pAdjustMethod %in% c(
        "bonferroni",
        "holm",
        "hochberg",
        "hommel",
        "BH")) {
stop("Invalid p-adjusted method. Please choose from 'bonferroni', 'holm',
'hochberg', 'hommel', or'fdr or BH'.")
    }

    # adjust p-value
    result[, "adjusted_pvalue"] <- p.adjust(result[, "pvalue"], method = "BH")
    result <- result[order(result[, "adjusted_pvalue"]), ]
    rownames(result) <- NULL
    result <- result[result$adjusted_pvalue <= pAdjustThreshold, ]

    # abs_b
    result <- result %>%
        mutate(abs_b = abs(result[, "b"]))

    result <- result[result$abs_b >= logfcThreshold, ]

    result_all <- rbind(result_all, result)
    }
    return(result_all)
}
