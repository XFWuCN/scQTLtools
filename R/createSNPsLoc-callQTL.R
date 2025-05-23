#' Create SNP location dataframe.
#'
#' @param snpList A list of SNPs IDs.
#' @param snp_mart  An object of class Mart representing the BioMart SNP
#' database to connect to. If provided, this should be a Mart object obtained
#' by calling \code{useEnsembl()}, which allows specifying a mirror in case of
#' connection issues. If \code{NULL}, the function will create and use a Mart
#' object pointing to the Ensembl SNP BioMart, using the specified
#' \code{snpDataset} and a default mirror.
#' @param snpDataset  A character string specifying the SNP dataset to use from
#' Ensembl. Default is \code{hsapiens_snp} for human SNPs.
#' @importFrom biomaRt useEnsembl getBM
#' @return A data frame containing the SNP genomic locations.
#' @export
#' @examples
#' snpList <- c('rs546', 'rs549', 'rs568', 'rs665', 'rs672')
#' library(biomaRt)
#' snp_mart <- useEnsembl(biomart = "snps",
#'                         dataset = "hsapiens_snp",
#'                         mirror = 'asia')
#' snp_loc <- createSNPsLoc(snpList = snpList,
#'                          snp_mart = snp_mart)
createSNPsLoc <- function(snpList,
                        snp_mart = NULL,
                        snpDataset = 'hsapiens_snp') {
    if(is.null(snp_mart)) {
        snp_mart <- useEnsembl(biomart = "snps", dataset = snpDataset)
    }
    stopifnot(is(snp_mart, 'Mart'))

    snps_loc <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                        filters = "snp_filter",
                        values = snpList,
                        mart = snp_mart)
    colnames(snps_loc)[which(colnames(snps_loc) == "chrom_start")] <-
        "position"
    rownames(snps_loc) <- snps_loc[, 1]
    return(snps_loc)
}
