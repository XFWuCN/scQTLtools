#' Check if the SNP ids in the input genotype matrix are valid.
#'
#' @param snpList a list of SNPs id.
#' @param snpDataset SNP dataset chosen from ENSEMBL
#' @return SNP location dataframe
#' @export
#' @examples
#' data(testSNP2)
#' snpList <- rownames(testSNP2)
#' snpDataset <- 'hsapiens_snp'
#' snps_loc <- checkSNPList(
#'   snpList = snpList,
#'   snpDataset = snpDataset
#'   )
checkSNPList <- function(snpList, snpDataset) {
    if (grepl("^rs", snpList[[1]][1])) {
        createSNPsLoc(snpList, snpDataset)
    } else if (grepl("\\d+:\\d+", snpList[[1]][1])) {
        snps_df <- data.frame(refsnp_id = character(),
                            chr_name = character(),
                            position = numeric(),
                            stringsAsFactors = FALSE)
        snps_loc <- snps_df
        for (i in seq_len(length(snpList))) {
            snp_parts <- strsplit(snpList[i], ":")[[1]]
            snps_loc <- rbind(snps_loc, data.frame(refsnp_id = snpList[i],
                chr_name = snp_parts[1], position = as.numeric(snp_parts[2])))
        }
        return(snps_loc)
    } else {
        return(NULL)
    }
}
