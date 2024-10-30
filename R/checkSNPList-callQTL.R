#' Check if the SNP ids in the input genotype matrix are valid.
#'
#' @param snpList a list of SNPs id.
#' @param snpDataset SNP dataset chosen from ENSEMBL
#' @param snpBiomart the name of snp mart.
#' @return SNP location dataframe
#' @export
#' @examples
#' data(testSNP2)
#' snpList <- rownames(testSNP2)
#' snpDataset <- 'hsapiens_snp'
#' snps_loc <- checkSNPList(snpList = snpList,
#'                         snpDataset = snpDataset)
checkSNPList <- function(snpList, snpDataset, snpBiomart) {
    if (grepl("^rs", snpList[[1]][1])) {
        createSNPsLoc(snpList, snpDataset, snpBiomart)
    } else if (grepl("\\d+:\\d+", snpList[[1]][1])) {
        snps_df <- data.frame(refsnp_id = character(),
                            chr_name = character(),
                            position = numeric(),
                            stringsAsFactors = FALSE)
        snps_loc_list <- lapply(snpList, function(snp) {
            snp_parts <- strsplit(snp, ":")[[1]]
            data.frame(refsnp_id = snp,
                        chr_name = snp_parts[1],
                        position = as.numeric(snp_parts[2]))
        })
        snps_loc1 <- do.call(rbind, snps_loc_list)
        return(snps_loc1)
    } else {
        return(NULL)
    }
}
