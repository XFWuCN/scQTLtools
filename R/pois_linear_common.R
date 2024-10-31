# @rdname poissonModel-utils_linearModel-utils_internals
process_snp_matrix <- function(snpMatrix_split, biClassify) {
    if (biClassify) {
        snpMatrix_split[snpMatrix_split == 3] <- 2
    } else {
    snpMatrix_split <- as.data.frame(snpMatrix_split) %>%
        mutate_all(~ifelse(. == 2, 3, ifelse(. == 3, 2, .)))
    snpMatrix_split <- as.matrix(snpMatrix_split)
    }
    return(snpMatrix_split)
}

# @rdname poissonModel-utils_linearModel-utils_internals
compute_gene_results_common <- function(snpID, geneIDs,
                    expressionMatrix_split, snpMatrix_split, group, model) {
    gene_results <- lapply(geneIDs, function(geneID) {
        snp_mat <- process_matrix(snpID, snpMatrix_split, "snp_mat")
        gene_mat <- process_matrix(geneID, expressionMatrix_split, "gene_mat")
        combined_df <- merge(snp_mat, gene_mat, by = "cells")
        combined_df <- subset(combined_df, snp_mat != 0)

    if (model == "linear") {
        lmodel <- lm(gene_mat ~ snp_mat, data = combined_df)
        if (length(summary(lmodel)$coefficients[, "Pr(>|t|)"]) >= 2) {
            return(data.frame(SNPid = snpID, group = group, Geneid = geneID,
                    pvalue = summary(lmodel)$coefficients[2, "Pr(>|t|)"],
                    b = summary(lmodel)$coefficients[2, "Estimate"]))
        }
    } else if (model == "poisson") {
        lmodel <- glm(combined_df$gene_mat ~ combined_df$snp_mat,
                        family = poisson())
        if (length(summary(lmodel)$coefficients[, "Pr(>|z|)"]) >= 2) {
            return(data.frame(SNPid = snpID, group = group, Geneid = geneID,
                pvalue = summary(lmodel)$coefficients[2, "Pr(>|z|)"],
                b = summary(lmodel)$coefficients[2, "Estimate"]))
            }
    }
    return(NULL)
    })
    return(do.call(rbind, gene_results))
}

# @rdname poissonModel-utils_linearModel-utils_internals
process_group_common <- function(eQTLObject, geneIDs, snpIDs, biClassify,
                                group, model) {
    split_cells <- rownames(load_group_info(eQTLObject)
                            )[load_group_info(eQTLObject)[["group"]] == group]
    expressionMatrix_split <- round(get_filter_data(eQTLObject)[["expMat"]] *
                            1000)[, split_cells, drop = FALSE]
    snpMatrix_split <- get_filter_data(eQTLObject)[["snpMat"]][, split_cells,
                                                            drop = FALSE]

    snpMatrix_split <- process_snp_matrix(snpMatrix_split, biClassify)
    pb <- initialize_progress_bar(length(snpIDs), group)

    result_list <- lapply(snpIDs, function(snpID) {
        result <- compute_gene_results_common(snpID, geneIDs,
                        expressionMatrix_split, snpMatrix_split, group, model)
        pb$tick()
        return(result)
    })
    return(do.call(rbind, result_list))
}
