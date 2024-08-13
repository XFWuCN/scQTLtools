#' Zinb model fitting the gene expression matrix.
#'
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param geneIDs Matching genes can be used to fit data.
#' @param snpIDs Matching SNPs can be used to fit data.
#' @param biClassify The user chooses whether to convert the counting method of
#' the snpMatrix to 0/1/2, TRUE indicates conversion, and FALSE indicates no
#' conversion, default is no conversion.
#' @param pAdjustMethod Methods for p-value adjusting, one of "bonferroni",
#' "holm", "hochberg", "hommel" or "BH". The default option is "bonferroni".
#' @param pAdjustThreshold  Only SNP-Gene pairs with adjusted p-values meeting
#' the threshold will be displayed. Default by 0.05.
#'
#' @importFrom stats pchisq plogis
#' @importFrom methods is
#' @importFrom VGAM dzinegbin
#' @return Dataframe that contains gene-SNP pairs' information.
#' @export
#' @examples
#' data(testEQTL)
#' Gene <- rownames(slot(testEQTL, "filterData")$normExpMat)
#' SNP <- rownames(slot(testEQTL, "filterData")$snpMat)
#' zinbResult <- zinbModel(
#'   eQTLObject = testEQTL,
#'   geneIDs = Gene,
#'   snpIDs = SNP,
#'   biClassify = FALSE,
#'   pAdjustMethod = "bonferroni",
#'   pAdjustThreshold = 0.05)
zinbModel <- function(
    eQTLObject,
    geneIDs,
    snpIDs,
    biClassify = FALSE,
    pAdjustMethod = "bonferroni",
    pAdjustThreshold = 1e-5) {
    if (length(eQTLObject@filterData) == 0) {
        stop("Please filter the data first.")
    }

    expressionMatrix <- round(eQTLObject@filterData$expMat * 1000)
    snpMatrix <- eQTLObject@filterData$snpMat
    unique_group <- unique(eQTLObject@groupBy$group)

    result_all <- data.frame()

    message("Start calling eQTL")
    for (j in unique_group) {
        split_cells <- rownames(eQTLObject@groupBy)[
            eQTLObject@groupBy$group == j]
        expressionMatrix_split <- expressionMatrix[, split_cells]
        snpMatrix_split <- snpMatrix[, split_cells]

    if (biClassify == TRUE) {
        snpMatrix_split[snpMatrix_split == 3] <- 2

        eQTLcalling <- function(i) {
            if (i %% 100 == 0) {
                gc()
            }

        snpid <- snpIDs[i]
        ref_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid, ] == 1]
        alt_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid, ] == 2]

        if (length(ref_cells) > 0 && length(alt_cells) > 0) {
            genes <- geneIDs
            gene.cnt <- 0

            results_snp <- data.frame(
                SNPid = character(),
                group = character(),
                Geneid = character(),
                sample_size_1 = integer(),
                sample_size_2 = integer(),
                theta_1 = double(),
                theta_2 = double(),
                mu_1 = double(),
                mu_2 = double(),
                size_1 = double(),
                size_2 = double(),
                prob_1 = double(),
                prob_2 = double(),
                total_mean_1 = double(),
                total_mean_2 = double(),
                foldChange = double(),
                chi = double(),
                pvalue = double(),
                adjusted_pvalue = double(),
                Remark = character(),
                stringsAsFactors = FALSE)

        for (gene in genes) {
            gene.cnt <- gene.cnt + 1
            counts_1 <- unlist(expressionMatrix_split[gene, ref_cells])
            counts_2 <- unlist(expressionMatrix_split[gene, alt_cells])
            results_gene <- data.frame(
                group = j,
                SNPid = snpid,
                Geneid = gene,
                sample_size_1 = length(counts_1),
                sample_size_2 = length(counts_2),
                theta_1 = NA,
                theta_2 = NA,
                mu_1 = NA,
                mu_2 = NA,
                size_1 = NA,
                size_2 = NA,
                prob_1 = NA,
                prob_2 = NA,
                total_mean_1 = NA,
                total_mean_2 = NA,
                foldChange = NA,
                chi = NA,
                pvalue = NA,
                adjusted_pvalue = NA,
                Remark = NA,
                stringsAsFactors = FALSE)

            totalMean_1 <- mean(counts_1)
            totalMean_2 <- mean(counts_2)
            foldChange <- totalMean_1 / totalMean_2

            results_gene[1, "total_mean_1"] <- totalMean_1
            results_gene[1, "total_mean_2"] <- totalMean_2
            results_gene[1, "foldChange"] <- foldChange


            params_1 <- buildZINB(counts_1)
            theta_1 <- params_1[["theta"]]
            mu_1 <- params_1[["mu"]]
            size_1 <- params_1[["size"]]
            prob_1 <- params_1[["prob"]]

            params_2 <- buildZINB(counts_2)
            theta_2 <- params_2[["theta"]]
            mu_2 <- params_2[["mu"]]
            size_2 <- params_2[["size"]]
            prob_2 <- params_2[["prob"]]

            params_combined <- buildZINB(c(counts_1, counts_2))
            theta_res <- params_combined[["theta"]]
            mu_res <- params_combined[["mu"]]
            size_res <- params_combined[["size"]]
            prob_res <- params_combined[["prob"]]

            logL <- function(counts_1,
                            theta_1,
                            size_1,
                            prob_1,
                            counts_2,
                            theta_2,
                            size_2,
                            prob_2) {
                logL_1 <- sum(
                    dzinegbin(counts_1,
                            size = size_1,
                            prob = prob_1,
                            pstr0 = theta_1,
                            log = TRUE))
                logL_2 <- sum(
                    dzinegbin(counts_2,
                            size = size_2,
                            prob = prob_2,
                            pstr0 = theta_2,
                            log = TRUE))
            logL <- logL_1 + logL_2
                logL
            }
            logL_A <- logL(
                counts_1,
                theta_1,
                size_1,
                prob_1,
                counts_2,
                theta_2,
                size_2,
                prob_2)
            logL_B <- logL(
                counts_1,
                theta_res,
                size_res,
                prob_res,
                counts_2,
                theta_res,
                size_res,
                prob_res)
            chi <- logL_A - logL_B
            pvalue <- 1 - pchisq(2 * chi, df = 3)

            results_gene[1, "theta_1"] <- theta_1
            results_gene[1, "theta_2"] <- theta_2
            results_gene[1, "mu_1"] <- mu_1
            results_gene[1, "mu_2"] <- mu_2
            results_gene[1, "size_1"] <- size_1
            results_gene[1, "size_2"] <- size_2
            results_gene[1, "prob_1"] <- prob_1
            results_gene[1, "prob_2"] <- prob_2
            results_gene[1, "chi"] <- chi
            results_gene[1, "pvalue"] <- pvalue

            results_snp <- rbind(results_snp, results_gene)
        }

            # return
            return(results_snp)
        } else {
            results_snp <- data.frame()
        }
    }


    # final result
    result <- data.frame(
        SNPid = character(),
        group = character(),
        Geneid = character(),
        sample_size_1 = integer(),
        sample_size_2 = integer(),
        theta_1 = double(),
        theta_2 = double(),
        mu_1 = double(),
        mu_2 = double(),
        size_1 = double(),
        size_2 = double(),
        prob_1 = double(),
        prob_2 = double(),
        total_mean_1 = double(),
        total_mean_2 = double(),
        foldChange = double(),
        chi = double(),
        pvalue = double(),
        adjusted_pvalue = double(),
        Remark = character(),
        stringsAsFactors = FALSE)

        message(j, ":")
        message("0%   10   20   30   40   50   60   70   80   90   100%")
        message("[----|----|----|----|----|----|----|----|----|----|")
        pb <- progress_bar$new(
            total = length(snpIDs),
            format = "[:bar]",
            clear = FALSE,
            width = 51)
    for (i in seq_len(length(snpIDs))) {
        result <- rbind(result, eQTLcalling(i))
        pb$tick()
    }
    colnames(result)[colnames(result) == "sample_size_1"] <- "sample_size_Ref"
    colnames(result)[colnames(result) == "sample_size_2"] <- "sample_size_Alt"
    colnames(result)[colnames(result) == "theta_1"] <- "theta_Ref"
    colnames(result)[colnames(result) == "theta_2"] <- "theta_Alt"
    colnames(result)[colnames(result) == "mu_1"] <- "mu_Ref"
    colnames(result)[colnames(result) == "mu_2"] <- "mu_Alt"
    colnames(result)[colnames(result) == "size_1"] <- "size_Ref"
    colnames(result)[colnames(result) == "size_2"] <- "size_Alt"
    colnames(result)[colnames(result) == "prob_1"] <- "prob_Ref"
    colnames(result)[colnames(result) == "prob_2"] <- "prob_Alt"
    colnames(result)[colnames(result) == "total_mean_1"] <- "total_mean_Ref"
    colnames(result)[colnames(result) == "total_mean_2"] <- "total_mean_Alt"
    } else if (biClassify == FALSE) {
        eQTLcalling <- function(i) {
            if (i %% 100 == 0) {
                gc()
            }

        snpid <- snpIDs[i]
        AA_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid, ] == 1]
        Aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid, ] == 3]
        aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid, ] == 2]

        if (length(AA_cells) > 0 &&
            length(Aa_cells) > 0 &&
            length(aa_cells) > 0) {
            genes <- geneIDs
            gene.cnt <- 0

        results_snp <- data.frame(
            SNPid = character(),
            group = character(),
            Geneid = character(),
            sample_size_1 = integer(),
            sample_size_2 = integer(),
            sample_size_3 = integer(),
            theta_1 = double(),
            theta_2 = double(),
            theta_3 = double(),
            mu_1 = double(),
            mu_2 = double(),
            mu_3 = double(),
            size_1 = double(),
            size_2 = double(),
            size_3 = double(),
            prob_1 = double(),
            prob_2 = double(),
            prob_3 = double(),
            total_mean_1 = double(),
            total_mean_2 = double(),
            total_mean_3 = double(),
            chi = double(),
            pvalue = double(),
            adjusted_pvalue = double(),
            Remark = character(),
            stringsAsFactors = FALSE)

        for (gene in genes) {
            gene.cnt <- gene.cnt + 1
            counts_1 <- unlist(expressionMatrix_split[gene, AA_cells])
            counts_2 <- unlist(expressionMatrix_split[gene, Aa_cells])
            counts_3 <- unlist(expressionMatrix_split[gene, aa_cells])
            results_gene <- data.frame(
                group = j,
                SNPid = snpid,
                Geneid = gene,
                sample_size_1 = length(counts_1),
                sample_size_2 = length(counts_2),
                sample_size_3 = length(counts_3),
                theta_1 = NA,
                theta_2 = NA,
                theta_3 = NA,
                mu_1 = NA,
                mu_2 = NA,
                mu_3 = NA,
                size_1 = NA,
                size_2 = NA,
                size_3 = NA,
                prob_1 = NA,
                prob_2 = NA,
                prob_3 = NA,
                total_mean_1 = NA,
                total_mean_2 = NA,
                total_mean_3 = NA,
                chi = NA,
                pvalue = NA,
                adjusted_pvalue = NA,
                Remark = NA,
                stringsAsFactors = FALSE)

            # calculate mean
            totalMean_1 <- mean(counts_1)
            totalMean_2 <- mean(counts_2)
            totalMean_3 <- mean(counts_3)

            # filling data
            results_gene[1, "total_mean_1"] <- totalMean_1
            results_gene[1, "total_mean_2"] <- totalMean_2
            results_gene[1, "total_mean_3"] <- totalMean_3

            # estimate AA params
            params_1 <- buildZINB(counts_1)
            theta_1 <- params_1[["theta"]]
            mu_1 <- params_1[["mu"]]
            size_1 <- params_1[["size"]]
            prob_1 <- params_1[["prob"]]

            # estimate Aa params
            params_2 <- buildZINB(counts_2)
            theta_2 <- params_2[["theta"]]
            mu_2 <- params_2[["mu"]]
            size_2 <- params_2[["size"]]
            prob_2 <- params_2[["prob"]]

            # estimate aa params
            params_3 <- buildZINB(counts_3)
            theta_3 <- params_3[["theta"]]
            mu_3 <- params_3[["mu"]]
            size_3 <- params_3[["size"]]
            prob_3 <- params_3[["prob"]]

            # combine AA,Aa and aa
            params_combined <- buildZINB(c(counts_1,
                                            counts_2,
                                            counts_3))
            theta_res <- params_combined[["theta"]]
            mu_res <- params_combined[["mu"]]
            size_res <- params_combined[["size"]]
            prob_res <- params_combined[["prob"]]

            # calculate p-value
                logL <- function(counts_1,
                                theta_1,
                                size_1,
                                prob_1,
                                counts_2,
                                theta_2,
                                size_2,
                                prob_2,
                                counts_3,
                                theta_3,
                                size_3,
                                prob_3) {
            logL_1 <- sum(dzinegbin(counts_1, size = size_1, prob = prob_1,
                pstr0 = theta_1, log = TRUE))
            logL_2 <- sum(dzinegbin(counts_2, size = size_2, prob = prob_2,
                pstr0 = theta_2, log = TRUE))
            logL_3 <- sum(dzinegbin(counts_3, size = size_3, prob = prob_3,
                pstr0 = theta_3, log = TRUE))
            logL <- logL_1 + logL_2 + logL_3
            logL
            }

            logL_1 <- logL(
                counts_1, theta_1, size_1, prob_1, counts_2, theta_2,
                size_2, prob_2, counts_3, theta_3, size_3, prob_3)
            logL_2 <- logL(
                counts_1, theta_res, size_res, prob_res, counts_2,
                theta_res, size_res, prob_res, counts_3, theta_res,
                size_res, prob_res)
            chi <- logL_1 - logL_2

            pvalue <- try(1 - pchisq(2 * chi, df = 3))
            if (is(pvalue, "try-error")) {
                return(NA)
            }

            # filling data into data frame
            results_gene[1, "theta_1"] <- theta_1
            results_gene[1, "theta_2"] <- theta_2
            results_gene[1, "theta_3"] <- theta_3
            results_gene[1, "mu_1"] <- mu_1
            results_gene[1, "mu_2"] <- mu_2
            results_gene[1, "mu_3"] <- mu_3
            results_gene[1, "size_1"] <- size_1
            results_gene[1, "size_2"] <- size_2
            results_gene[1, "size_3"] <- size_3
            results_gene[1, "prob_1"] <- prob_1
            results_gene[1, "prob_2"] <- prob_2
            results_gene[1, "prob_3"] <- prob_3
            results_gene[1, "chi"] <- chi
            results_gene[1, "pvalue"] <- pvalue

            results_snp <- rbind(results_snp, results_gene)
            }

        # return
        return(results_snp)
        } else {
            results_snp <- data.frame()
        }
    }

    # final result
    result <- data.frame(
        group = character(),
        SNPid = character(),
        Geneid = character(),
        sample_size_1 = integer(),
        sample_size_2 = integer(),
        sample_size_3 = integer(),
        theta_1 = double(),
        theta_2 = double(),
        theta_3 = double(),
        mu_1 = double(),
        mu_2 = double(),
        mu_3 = double(),
        size_1 = double(),
        size_2 = double(),
        size_3 = double(),
        prob_1 = double(),
        prob_2 = double(),
        prob_3 = double(),
        total_mean_1 = double(),
        total_mean_2 = double(),
        total_mean_3 = double(),
        chi = double(),
        pvalue = double(),
        adjusted_pvalue = double(),
        Remark = character(),
        stringsAsFactors = FALSE)

        message(j, ":")
        message("0%   10   20   30   40   50   60   70   80   90   100%")
        message("[----|----|----|----|----|----|----|----|----|----|")
        pb <- progress_bar$new(
            total = length(snpIDs),
            format = "[:bar]",
            clear = FALSE,
            width = 51)
    for (i in seq_len(length(snpIDs))) {
        result <- rbind(result, eQTLcalling(i))
        pb$tick()
        }

        message("finished!")

    # change column names
    colnames(result)[colnames(result) == "sample_size_1"] <- "sample_size_AA"
    colnames(result)[colnames(result) == "sample_size_2"] <- "sample_size_Aa"
    colnames(result)[colnames(result) == "sample_size_3"] <- "sample_size_aa"
    colnames(result)[colnames(result) == "theta_1"] <- "theta_AA"
    colnames(result)[colnames(result) == "theta_2"] <- "theta_Aa"
    colnames(result)[colnames(result) == "theta_3"] <- "theta_aa"
    colnames(result)[colnames(result) == "mu_1"] <- "mu_AA"
    colnames(result)[colnames(result) == "mu_2"] <- "mu_Aa"
    colnames(result)[colnames(result) == "mu_3"] <- "mu_aa"
    colnames(result)[colnames(result) == "size_1"] <- "size_AA"
    colnames(result)[colnames(result) == "size_2"] <- "size_Aa"
    colnames(result)[colnames(result) == "size_3"] <- "size_aa"
    colnames(result)[colnames(result) == "prob_1"] <- "prob_AA"
    colnames(result)[colnames(result) == "prob_2"] <- "prob_Aa"
    colnames(result)[colnames(result) == "prob_3"] <- "prob_aa"
    colnames(result)[colnames(result) == "total_mean_1"] <- "total_mean_AA"
    colnames(result)[colnames(result) == "total_mean_2"] <- "total_mean_Aa"
    colnames(result)[colnames(result) == "total_mean_3"] <- "total_mean_aa"
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
    result[, "adjusted_pvalue"] <- p.adjust(result[, "pvalue"],
                                            method = pAdjustMethod)
    result <- result[order(result[, "adjusted_pvalue"]), ]
    rownames(result) <- NULL
    result <- result[result$adjusted_pvalue <= pAdjustThreshold, ]
    result_all <- rbind(result_all, result)
    }
    return(result_all)
}
