#' Create gene location dataframe.
#'
#' @param geneList A gene id or a list of genes id.
#' @param geneDataset Gene dataset chosen from the biomart.
#' @param OrgDb OrgDb name:"org.Hs.eg.db", "org.Mm.eg.db", "org.Ce.eg.db",
#' "org.At.tair.db".
#' @param geneBiomart the name of gene mart.
#' @param host Host to connect to.
#' @param gene_attributes the outputs of a biomaRt query.
#' @param filters what we use as inputs for a biomaRt query.
#' @importFrom biomaRt getBM useEnsembl
#' @importFrom AnnotationDbi mapIds
#' @return data.frame
#' @export
#' @examples
#' data(testGene)
#' geneList <- rownames(testGene)
#' geneDataset <- 'hsapiens_gene_ensembl'
#' library(GOSemSim)
#' OrgDb <- load_OrgDb("org.Hs.eg.db")
#' gene_loc <- createGeneLoc(geneList = geneList,
#'                           geneDataset = geneDataset,
#'                           OrgDb = OrgDb)
createGeneLoc <- function(geneList,
                            geneDataset,
                            OrgDb,
                            geneBiomart = "ENSEMBL_MART_ENSEMBL",
                            host = "https://www.ensembl.org",
                            gene_attributes = c("external_gene_name",
                                                "chromosome_name",
                                                "start_position",
                                                "end_position"),
                            filters = "external_gene_name") {
    gene_mart <- useMart(biomart = geneBiomart,
                        dataset = geneDataset,
                        host = host)
    geneList <- unique(geneList)
    if (grepl("^ENSG", geneList[[1]][1])) {
        gene_attributes <- c("ensembl_gene_id",
                            "chromosome_name",
                            "start_position",
                            "end_position")
        filters <- "ensembl_gene_id"
    } else {
        gene_attributes <- gene_attributes
        filters <- filters
    }
    gene_loc <- getBM(attributes = gene_attributes,
                        filters = filters,
                        values = geneList,
                        mart = gene_mart)
    gene_loc <- gene_loc[grepl("^[0-9]+$",
                        as.character(gene_loc$chromosome_name)), ]
    rownames(gene_loc) <- NULL
    return(gene_loc)
}
