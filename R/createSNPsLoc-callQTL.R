#' Create SNP location dataframe.
#'
#' @param snpList a list of SNPs id.
#' @param snpDataset SNP dataset chosen from ENSEMBL
#' @param snpBiomart the name of snp mart.
#' @importFrom biomaRt useMart getBM
#' @return data.frame
#' @export
#' @examples
#' snpList <- c('rs546', 'rs549', 'rs568', 'rs665', 'rs672')
#' snpDataset <- 'hsapiens_snp'
#' snpBiomart <- "ENSEMBL_MART_SNP"
#' snp_loc <- createSNPsLoc(snpList = snpList,
#'                          snpDataset = snpDataset,
#'                          snpBiomart = snpBiomart)
createSNPsLoc <- function(snpList,
                        snpDataset,
                        snpBiomart = "ENSEMBL_MART_SNP") {
    snp_mart <- useMart(biomart = snpBiomart,
                        dataset = snpDataset)
    snps_loc <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                        filters = "snp_filter",
                        values = snpList,
                        mart = snp_mart)
    colnames(snps_loc)[which(colnames(snps_loc) == "chrom_start")] <-
        "position"
    rownames(snps_loc) <- snps_loc[, 1]
    return(snps_loc)
}
