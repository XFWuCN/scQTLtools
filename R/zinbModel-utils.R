#' Zinb model fitting the gene expression matrix.
#'
#' @param eQTLObject An S4 object of class \code{eQTLObject}.
#' @param geneIDs Character vector of gene IDs to include in the model fitting.
#' @param snpIDs Character vector of SNP IDs to include in the model fitting.
#' @param biClassify Logical; whether to convert genotype encoding in snpMatrix
#' to 0, 1, and 2. \code{TRUE} indicates conversion; \code{FALSE} indicates no
#' conversion (default).
#' @param pAdjustThreshold  Only gene-SNP pairs with adjusted p-values below
#' the threshold will be retained. Default is 0.05.
#' @param pAdjustMethod  Method used for multiple testing correction. One of
#' \code{"bonferroni"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, or
#' \code{"BH"}. Default is \code{"bonferroni"}.
#'
#' @importFrom stats pchisq plogis
#' @importFrom methods is
#' @importFrom VGAM dzinegbin
#' @return A data frame of gene–SNP pairs that pass the filtering criteria,
#' including p-values, adjusted p-values, and group labels.
#' @export
#' @examples
#' data(testEQTL)
#' Gene <- rownames(slot(testEQTL, 'filterData')$expMat)
#' SNP <- rownames(slot(testEQTL, 'filterData')$snpMat)
#' zinbResult <- zinbModel(
#'   eQTLObject = testEQTL,
#'   geneIDs = Gene,
#'   snpIDs = SNP,
#'   biClassify = FALSE,
#'   pAdjustMethod = 'bonferroni',
#'   pAdjustThreshold = 0.05)
zinbModel <- function(eQTLObject, geneIDs, snpIDs, biClassify = FALSE,
    pAdjustMethod = "bonferroni", pAdjustThreshold = 0.05) {
    expressionMatrix <- round(get_filter_data(eQTLObject)[["expMat"]] * 1000)
    snpMatrix <- get_filter_data(eQTLObject)[["snpMat"]]
    unique_group <- unique(load_group_info(eQTLObject)[["group"]])
    message("Start the sc-eQTLs calling.")
    result_all <- do.call(rbind, lapply(unique_group, function(j) {
        process_group(eQTLObject, j, expressionMatrix, snpMatrix, geneIDs,
            snpIDs, biClassify)
    }))
    result_all <- adjust_pvalues(result_all, pAdjustMethod, pAdjustThreshold)
    return(result_all)
}

# @rdname zinbModel-utils_internals
process_group <- function(eQTLObject, group, expressionMatrix, snpMatrix,
                        geneIDs, snpIDs, biClassify) {
    split_cells <- rownames(load_group_info(eQTLObject))[
        load_group_info(eQTLObject)[["group"]] == group]
    expressionMatrix_split <- expressionMatrix[, split_cells, drop = FALSE]
    snpMatrix_split <- snpMatrix[, split_cells, drop = FALSE]
    if (biClassify) {
        snpMatrix_split[snpMatrix_split == 3] <- 2
    }
    pb <- initialize_progress_bar(length(snpIDs), group)
    results_list <- lapply(seq_along(snpIDs), function(i) {
        result <- eQTLcalling(i, snpIDs, expressionMatrix_split,
                            snpMatrix_split, geneIDs, biClassify)
        pb$tick()
        return(result)
    })
    result <- do.call(rbind, results_list)
    result <- rename_columns(result, biClassify)
    result$group <- group
    return(result)
}
# @rdname zinbModel-utils_internals
eQTLcalling <- function(index, snpIDs, expressionMatrix_split, snpMatrix_split,
    geneIDs, biClassify) {
    snpid <- snpIDs[index]
    cell_groups <- get_cell_groups(snpMatrix_split, snpid, biClassify)
    if (length(cell_groups) < 2)
        return(data.frame())
    counts_list <- lapply(geneIDs, function(gene) {
        counts <- lapply(cell_groups,
            function(cells) get_counts(expressionMatrix_split, gene, cells))
        calculate_results(snpid, gene, counts)
    })
    return(do.call(rbind, counts_list))
}
# @rdname zinbModel-utils_internals
calculate_results <- function(snpid, gene, counts) {
    counts_1 <- counts[[1]]
    counts_2 <- counts[[2]]
    if (length(counts_1) == 0 || length(counts_2) == 0)
        return(data.frame())
    results_gene <- initialize_results_dataframe(snpid, gene, counts)
    params_1 <- buildZINB(counts_1)
    params_2 <- buildZINB(counts_2)
    params_combined <- buildZINB(c(counts_1, counts_2))
    results_gene <- extract_parameters(results_gene, params_1, params_2,
        params_combined)
    chi <- calculate_chi(counts_1, counts_2, params_1, params_2,
        params_combined)
    pvalue <- 1 - pchisq(2 * chi, df = 3)
    results_gene[1, "pvalue"] <- pvalue
    results_gene[1, "chi"] <- chi
    return(results_gene)
}
# @rdname zinbModel-utils_internals
initialize_results_dataframe <- function(snpid, gene, counts) {
    return(data.frame(SNPid = snpid, Geneid = gene,
        sample_size_1 = length(counts[[1]]),
        sample_size_2 = length(counts[[2]]),
        theta_1 = NA, theta_2 = NA, mu_1 = NA,
        mu_2 = NA, size_1 = NA, size_2 = NA,
        prob_1 = NA, prob_2 = NA, total_mean_1 = mean(counts[[1]],
            na.rm = TRUE), total_mean_2 = mean(counts[[2]],
            na.rm = TRUE), foldChange = mean(counts[[1]],
            na.rm = TRUE)/mean(counts[[2]], na.rm = TRUE),
        chi = NA, pvalue = NA, adjusted_pvalue = NA,
        Remark = NA, stringsAsFactors = FALSE))
}
# @rdname zinbModel-utils_internals
extract_parameters <- function(results_gene, params_1, params_2,
    params_combined) {
    results_gene[1, "theta_1"] <- params_1[["theta"]]
    results_gene[1, "theta_2"] <- params_2[["theta"]]
    results_gene[1, "mu_1"] <- params_1[["mu"]]
    results_gene[1, "mu_2"] <- params_2[["mu"]]
    results_gene[1, "size_1"] <- params_1[["size"]]
    results_gene[1, "size_2"] <- params_2[["size"]]
    results_gene[1, "prob_1"] <- params_1[["prob"]]
    results_gene[1, "prob_2"] <- params_2[["prob"]]
    return(results_gene)
}
# @rdname zinbModel-utils_internals
calculate_chi <- function(counts_1, counts_2, params_1, params_2,
    params_combined) {
    logL_A <- logL(counts_1, params_1, counts_2, params_2)
    logL_B <- logL(counts_1, params_combined, counts_2, params_combined)
    return(logL_A - logL_B)
}
# @rdname zinbModel-utils_internals
logL <- function(counts_1, params_1, counts_2, params_2) {
    logL_1 <- sum(dzinegbin(counts_1, size = params_1[["size"]],
        prob = params_1[["prob"]], pstr0 = params_1[["theta"]], log = TRUE))
    logL_2 <- sum(dzinegbin(counts_2, size = params_2[["size"]],
        prob = params_2[["prob"]], pstr0 = params_2[["theta"]], log = TRUE))
    return(logL_1 + logL_2)
}
# @rdname zinbModel-utils_internals
rename_columns <- function(result, biClassify) {
    if (biClassify) {
        colnames(result)[colnames(result) ==
            "sample_size_1"] <- "sample_size_Ref"
        colnames(result)[colnames(result) ==
            "sample_size_2"] <- "sample_size_Alt"
        colnames(result)[colnames(result) ==
            "theta_1"] <- "theta_Ref"
        colnames(result)[colnames(result) ==
            "theta_2"] <- "theta_Alt"
    } else {
        colnames(result)[colnames(result) ==
            "sample_size_1"] <- "sample_size_AA"
        colnames(result)[colnames(result) ==
            "sample_size_2"] <- "sample_size_Aa"
        colnames(result)[colnames(result) ==
            "sample_size_3"] <- "sample_size_aa"
        colnames(result)[colnames(result) ==
            "theta_1"] <- "theta_AA"
        colnames(result)[colnames(result) ==
            "theta_2"] <- "theta_Aa"
        colnames(result)[colnames(result) ==
            "theta_3"] <- "theta_aa"
    }
    return(result)
}
