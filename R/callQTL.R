#' callQTL: Uncover single-cell eQTLs exclusively using scRNA-seq data.
#' A function designed to identify eQTLs from scRNA-seq data.
#' @param useModel Model for fitting dataframe, one of 'possion', 'zinb', or
#' 'linear'.
#' @param pAdjustThreshold Only SNP-Gene pairs with adjusted p-values meeting
#' the threshold will be displayed. Default by 0.05.
#' @param pAdjustMethod Methods for p-value adjusting, one of 'bonferroni',
#' 'holm', 'hochberg', 'hommel' or 'BH'. Default by 'bonferroni'.
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param gene_ids A gene ID or a list of gene IDS.
#' @param downstream Being used to match SNPs within a base range defined by
#' the start position of genes.
#' @param upstream Being used to match SNPs within a base range defined by the
#' end position of genes.
#' @param logfcThreshold Represents the minimum beta threshold for fitting
#' SNP-Gene pairs.
#' @importFrom Matrix Matrix
#' @importFrom stringr str_split
#' @importFrom dplyr mutate_all mutate
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import magrittr
#' @return A dataframe, each row describes eQTL discovering result of a
#' SNP-Gene pair.
#' @export
#' @examples
#' data(testEQTL)
#' eqtl <- callQTL(
#'   eQTLObject = testEQTL,
#'   gene_ids = NULL,
#'   downstream = NULL,
#'   upstream = NULL,
#'   pAdjustMethod = 'bonferroni',
#'   useModel = 'linear',
#'   pAdjustThreshold = 0.05,
#'   logfcThreshold = 0.025
#' )
callQTL <- function(eQTLObject,
                    gene_ids = NULL,
                    downstream = NULL,
                    upstream = NULL,
                    pAdjustMethod = "bonferroni",
                    useModel = "zinb",
                    pAdjustThreshold = 0.05,
                    logfcThreshold = 0.1) {
    options(warn = -1)

    if (length(eQTLObject@filterData) == 0) {
        stop("Please filter the data first.")
    } else {
    expressionMatrix <- eQTLObject@filterData$expMat
    snpMatrix <- eQTLObject@filterData$snpMat
    }

    biClassify <- eQTLObject@biClassify
    species <- eQTLObject@species

    if (is.null(gene_ids) && is.null(upstream) && is.null(downstream)) {
    NULL
    } else {
    if (!is.null(species) && species != "") {
        if (species == "human") {
            snpDataset <- "hsapiens_snp"
            geneDataset <- "hsapiens_gene_ensembl"
            OrgDb <- org.Hs.eg.db
    } else if (species == "mouse") {
            snpDataset <- "mmusculus_snp"
            geneDataset <- "mmusculus_gene_ensembl"
            OrgDb <- org.Mm.eg.db
        } else {
            stop("Please enter 'human' or 'mouse'.")
        }
    } else {
        stop("The 'species' variable is NULL or empty.")
        }
    }

    snpList <- rownames(snpMatrix)
    geneList <- rownames(expressionMatrix)

    if (is.null(gene_ids) && is.null(upstream) && is.null(downstream)) {
        matched_gene <- geneList
        matched_snps <- snpList
    } else if (!is.null(gene_ids) && is.null(upstream) && is.null(downstream)) {
        matched_snps <- snpList
    if (all(gene_ids %in% geneList)) {
        matched_gene <- gene_ids
    } else {
    stop("The input gene_ids contain non-existent gene IDs.Please re-enter.")
    }
    } else if (is.null(gene_ids) &&
                !is.null(upstream) &&
                !is.null(downstream)) {
    if (downstream > 0) {
        stop("downstream should be negative.")
    }
    snps_loc <- checkSNPList(snpList, snpDataset)
    gene_loc <- createGeneLoc(geneList, geneDataset)

    matched_gene <- c()
    matched_snps <- c()

    for (i in seq_len(nrow(gene_loc))) {
        gene_start1 <- gene_loc$start_position[i] + downstream
        gene_end1 <- gene_loc$end_position[i] + upstream

        for (j in seq_len(nrow(snps_loc))) {
            snp_chr <- snps_loc$chr_name[j]
            snp_pos <- snps_loc$position[j]

        if (snp_chr == gene_loc$chromosome_name[i] &&
            snp_pos >= gene_start1 && snp_pos <= gene_end1) {
                matched_snps <- c(matched_snps, snps_loc$refsnp_id[j])
                matched_gene <- c(matched_gene, gene_loc[i, 1])
            }
        }
    }
    matched_snps <- unique(matched_snps)
    matched_gene <- unique(matched_gene)
    } else if (!is.null(gene_ids) &&
                !is.null(upstream) &&
                !is.null(downstream)) {
    if (downstream > 0) {
        stop("downstream should be negative.")
    }

    snps_loc <- checkSNPList(snpList)

    if (all(gene_ids %in% geneList)) {
        geneList <- gene_ids
    } else {
    stop("The input gene_ids contain non-existent gene IDs.Please re-enter.")
    }
    gene_loc <- createGeneLoc(geneList)

    matched_gene <- c()
    matched_snps <- c()

    for (i in seq_len(nrow(gene_loc))) {
        gene_start1 <- gene_loc$start_position[i] + downstream
        gene_end1 <- gene_loc$end_position[i] + upstream

        for (j in seq_len(nrow(snps_loc))) {
            snp_chr <- snps_loc$chr_name[j]
            snp_pos <- snps_loc$position[j]

        while (snp_chr == gene_loc$chromosome_name[i] &&
            snp_pos >= gene_start1 && snp_pos <= gene_end1) {
                matched_snps <- c(matched_snps, snps_loc$refsnp_id[j])
                matched_gene <- c(matched_gene, gene_loc[i, 1])
            }
        }
    }
    matched_snps <- unique(matched_snps)
    matched_gene <- unique(matched_gene)
    } else {
        stop("Please enter upstream and downstream simultaneously.")
        }


    if (useModel == "zinb") {
    result <- zinbModel(eQTLObject = eQTLObject,
                        geneIDs = matched_gene,
                        snpIDs = matched_snps,
                        biClassify = biClassify,
                        pAdjustMethod = pAdjustMethod,
                        pAdjustThreshold = pAdjustThreshold)
    } else if (useModel == "poisson") {
    result <- poissonModel(
        eQTLObject = eQTLObject,
        geneIDs = matched_gene,
        snpIDs = matched_snps,
        biClassify = biClassify,
        pAdjustMethod = pAdjustMethod,
        pAdjustThreshold = pAdjustThreshold,
        logfcThreshold = logfcThreshold
        )
    } else if (useModel == "linear") {
    result <- linearModel(
        eQTLObject = eQTLObject,
        geneIDs = matched_gene,
        snpIDs = matched_snps,
        biClassify = biClassify,
        pAdjustMethod = pAdjustMethod,
        pAdjustThreshold = pAdjustThreshold,
        logfcThreshold = logfcThreshold)
    } else {
        stop("Invalid model Please choose from 'zinb','poisson',or 'linear'.")
        }

    options(warn = 0)
    eQTLObject@useModel <- useModel
    eQTLObject@eQTLResult <- result
    return(eQTLObject)
    }
