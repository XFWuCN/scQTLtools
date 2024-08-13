#' Create gene location dataframe
#'
#' @param geneList a gene id or a list of genes id
#' @param geneDataset Gene dataset chosen from the biomart.
#' @importFrom biomaRt getBM useEnsembl
#' @importFrom AnnotationDbi mapIds
#' @return data.frame
#' @export
#' @examples
#' data(testGene)
#' geneList <- rownames(testGene)
#' geneDataset <- 'hsapiens_gene_ensembl'
#' gene_loc <- createGeneLoc(geneList = geneList,
#'                           geneDataset = geneDataset)
createGeneLoc <- function(geneList, geneDataset) {
    gene_mart <- useEnsembl(biomart = "ensembl",
                            dataset = geneDataset)
    geneList <- unique(geneList)

    if (grepl("^ENSG", geneList[[1]][1])) {
    gene_name <- "ensembl_gene_id"
    } else {
    gene_name <- "external_gene_name"
    ensembls <- mapIds(
        org.Hs.eg.db,
        keys = geneList,
        keytype = "SYMBOL",
        column = "ENSEMBL",
        multiVals = "first"
        )
    ensembls <- as.data.frame(ensembls)
    ensembls_id <- unique(ensembls$ensembls)
    }

    gene_attributes <- c(gene_name,
                        "chromosome_name",
                        "start_position",
                        "end_position")
    gene_loc <- getBM(attributes = gene_attributes,
                        filters = gene_name,
                        values = geneList,
                        mart = gene_mart)

    gene_loc <- gene_loc[grepl("^[0-9]+$",
                        as.character(gene_loc$chromosome_name)), ]
    rownames(gene_loc) <- NULL
    return(gene_loc)
}
