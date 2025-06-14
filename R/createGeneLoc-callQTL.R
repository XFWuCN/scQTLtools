#' Create a gene location data frame.
#'
#' @param geneList A gene ID or a list of genes IDs.
#' @param geneDataset  A character string specifying the gene dataset to use
#' from Ensembl. Default is \code{hsapiens_gene_ensembl} for human genes.
#' @param OrgDb  The name of the OrgDb package to use for gene annotation.
#' Supported values include \code{org.Hs.eg.db} and \code{org.Mm.eg.db}.
#' @param gene_mart  An object of class Mart representing the BioMart gene
#' database to connect to. If provided, this should be a Mart object obtained
#' by calling \code{useEnsembl()}, which allows specifying a mirror in case of
#' connection issues. If \code{NULL}, the function will create and use a Mart
#' object pointing to the Ensembl Gene BioMart, using the specified
#' \code{geneDataset} and a default mirror.
#' @importFrom biomaRt getBM useEnsembl
#' @return A \code{data.frame} containing gene location information.
#' @export
#' @examples
#' data(GeneData)
#' geneList <- rownames(GeneData)
#' library(GOSemSim)
#' library(biomaRt)
#' OrgDb <- load_OrgDb("org.Hs.eg.db")
#' gene_mart <- useEnsembl(biomart = "genes",
#'                         dataset = "hsapiens_gene_ensembl",
#'                         mirror = 'asia')
#' gene_loc <- createGeneLoc(geneList = geneList,
#'                           gene_mart = gene_mart,
#'                           OrgDb = OrgDb)
createGeneLoc <- function(geneList,
                        gene_mart = NULL,
                        geneDataset = "hsapiens_gene_ensembl",
                        OrgDb) {
    if(is.null(gene_mart)) {
        gene_mart <- useEnsembl(biomart = "genes", dataset = geneDataset)
    }
    stopifnot(is(gene_mart, 'Mart'))

    geneList <- unique(geneList)

    if (grepl("^ENSG", geneList[[1]][1])) {
        gene_name <- "ensembl_gene_id"
    } else {
        gene_name <- "external_gene_name"
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
