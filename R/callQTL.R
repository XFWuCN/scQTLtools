#' Identify single-cell eQTLs.
#'
#' This function detects eQTLs using scRNA-seq data and genotype data.
#' @param useModel  Model used for fitting. One of "poisson", "zinb", or
#' "linear."
#' @param pAdjustThreshold  Only SNP–gene pairs with adjusted p-values below
#' the threshold will be retained. Default is 0.05.
#' @param pAdjustMethod  Method used for multiple testing correction. One of
#' \code{"bonferroni"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, or
#' \code{"BH"}. Default is \code{"bonferroni"}.
#' @param eQTLObject  An S4 object of class \code{eQTLObject}.
#' @param gene_ids  A gene ID or a list of gene IDS.
#' @param downstream  Distance (in base pairs) downstream of the gene start
#' site to search for associated SNPs.
#' @param upstream  Distance (in base pairs) upstream of the gene end site to
#' search for associated SNPs.
#' @param gene_mart  A Mart object representing the BioMart gene database.
#' If \code{NULL}, the Ensembl Gene BioMart will be used.
#' @param snp_mart  A Mart object representing the BioMart SNP database.
#' If \code{NULL}, the Ensembl SNP BioMart will be used.
#' @param logfcThreshold  The minimum beta coefficient required to report a
#' SNP–gene pair as an eQTL.
#' @importFrom Matrix Matrix
#' @importFrom stringr str_split
#' @importFrom dplyr mutate_all mutate
#' @importFrom GOSemSim load_OrgDb
#' @importFrom stats na.omit
#' @import magrittr
#' @return A data frame in which each row corresponds to a detected SNP–gene
#' eQTL pair, including statistical and model fitting results.
#' @export
#' @examples
#' data(EQTL_obj)
#' library(biomaRt)
#' gene_mart <- useEnsembl(biomart = "genes",
#'                         dataset = "hsapiens_gene_ensembl",
#'                         mirror = 'asia')
#' snp_mart <- useEnsembl(biomart = "snps",
#'                         dataset = "hsapiens_snp",
#'                         mirror = 'asia')
#' eqtl <- callQTL(
#'   eQTLObject = EQTL_obj,
#'   gene_ids = NULL,
#'   downstream = NULL,
#'   upstream = NULL,
#'   gene_mart = gene_mart,
#'   snp_mart = snp_mart,
#'   pAdjustMethod = 'bonferroni',
#'   useModel = 'linear',
#'   pAdjustThreshold = 0.05,
#'   logfcThreshold = 0.025
#' )
callQTL <- function(eQTLObject,
                    gene_ids = NULL,
                    downstream = NULL,
                    upstream = NULL,
                    gene_mart = NULL,
                    snp_mart = NULL,
                    pAdjustMethod = "bonferroni",
                    useModel = "zinb",
                    pAdjustThreshold = 0.05,
                    logfcThreshold = 0.1) {
    options(warn = -1)
    if (length(get_filter_data(eQTLObject)) == 0) {
        stop("Please filter the data first.")
    } else {
        expressionMatrix <- get_filter_data(eQTLObject)[["expMat"]]
        snpMatrix <- get_filter_data(eQTLObject)[["snpMat"]]
    }
    biClassify <- load_biclassify_info(eQTLObject)
    species <- load_species_info(eQTLObject)
    if (is.null(species) || species == "") {
        stop("The 'species' variable is NULL or empty.")
    }
    if (species == "human") {
        snpDataset <- "hsapiens_snp"
        geneDataset <- "hsapiens_gene_ensembl"
        OrgDb <- "org.Hs.eg.db"
    } else if (species == "mouse") {
        snpDataset <- "mmusculus_snp"
        geneDataset <- "mmusculus_gene_ensembl"
        OrgDb <- "org.Mm.eg.db"
    } else {
        stop("Please enter 'human' or 'mouse'.")
    }
    if (is.null(gene_ids) && is.null(upstream) && is.null(downstream)) {
        NULL
    }
    snpList <- rownames(snpMatrix)
    geneList <- rownames(expressionMatrix)
    matchID <- match_gene_snp(gene_ids, upstream, downstream, snpList,
            geneList, gene_mart, snp_mart, geneDataset, snpDataset, OrgDb)
    matched_gene <- matchID[["matched_gene"]]
    matched_snps <- matchID[["matched_snps"]]
    result <- run_model_result(eQTLObject = eQTLObject, geneIDs = matched_gene,
    snpIDs = matched_snps, useModel = useModel, pAdjustMethod = pAdjustMethod,
    pAdjustThreshold = pAdjustThreshold, logfcThreshold = logfcThreshold)
    options(warn = 0)
    eQTLObject <- set_model_info(eQTLObject, useModel)
    eQTLObject <- set_result_info(eQTLObject, result)
    return(eQTLObject)
    }


# @rdname callQTL_internals
match_gene_snp <- function(gene_ids, upstream, downstream, snpList, geneList,
                            gene_mart = NULL, snp_mart = NULL, geneDataset,
                            snpDataset, OrgDb){
    if (is.null(gene_ids) && is.null(upstream) && is.null(downstream)) {
        matched_gene <- geneList
        matched_snps <- snpList
    } else if (!is.null(gene_ids) && is.null(upstream) &&
                is.null(downstream)) {
        matched_snps <- snpList
    if (all(gene_ids %in% geneList)) {
        matched_gene <- gene_ids
    } else {
    stop("The input gene_ids contain non-existent gene IDs. Please re-enter.")
    }
    } else if (is.null(gene_ids) && !is.null(upstream) &&
                !is.null(downstream)) {
    if (downstream > 0) {
        stop("downstream should be negative.")
    }

    snps_loc <- checkSNPList(snpList, snp_mart, snpDataset)
    gene_loc <- createGeneLoc(geneList, gene_mart, geneDataset, OrgDb)
    gene_ranges <- data.frame(
        gene_start = gene_loc$start_position + downstream,
        gene_end = gene_loc$end_position + upstream,
        chr = gene_loc$chromosome_name,
        gene_id = gene_loc[, 1]
    )
    matches <- lapply(seq_len(nrow(gene_ranges)), function(i) {
        snp_matches <- snps_loc[snps_loc$chr_name == gene_ranges$chr[i] &
                        snps_loc$position >= gene_ranges$gene_start[i] &
                        snps_loc$position <= gene_ranges$gene_end[i], ]
        if (nrow(snp_matches) > 0) {
            return(data.frame(snp_id = snp_matches$refsnp_id,
                                gene = gene_ranges$gene_id[i]))
        } else {
            return(NULL)
        }
    })
        matches_df <- do.call(rbind, matches)
        matches_df <- na.omit(matches_df)
        matched_snps <- unique(matches_df$snp_id)
        matched_gene <- unique(matches_df$gene)
    } else {
        stop("Please enter upstream and downstream simultaneously.")
    }
    return(list(matched_gene = matched_gene, matched_snps = matched_snps))
}


# @rdname callQTL_internals
run_model_result <- function(eQTLObject,
                            geneIDs,
                            snpIDs,
                            useModel,
                            pAdjustMethod,
                            pAdjustThreshold,
                            logfcThreshold){
    biClassify <- load_biclassify_info(eQTLObject)
    if (useModel == "zinb") {
        result <- zinbModel(
            eQTLObject = eQTLObject,
            geneIDs = geneIDs,
            snpIDs = snpIDs,
            biClassify = biClassify,
            pAdjustMethod = pAdjustMethod,
            pAdjustThreshold = pAdjustThreshold)
    } else if (useModel == "poisson") {
        result <- poissonModel(
            eQTLObject = eQTLObject,
            geneIDs = geneIDs,
            snpIDs = snpIDs,
            biClassify = biClassify,
            pAdjustMethod = pAdjustMethod,
            pAdjustThreshold = pAdjustThreshold,
            logfcThreshold = logfcThreshold)
    } else if (useModel == "linear") {
        result <- linearModel(
            eQTLObject = eQTLObject,
            geneIDs = geneIDs,
            snpIDs = snpIDs,
            biClassify = biClassify,
            pAdjustMethod = pAdjustMethod,
            pAdjustThreshold = pAdjustThreshold,
            logfcThreshold = logfcThreshold)
    } else {
        stop("Invalid model Please choose from 'zinb','poisson',or 'linear'.")
    }
    return(result)
}
