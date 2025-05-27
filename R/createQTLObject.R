# Defines a new class union named "AnyMatrixOrDataframe."
# This union allows objects of class "matrix", "dgCMatrix", or
# "data.frame" to be treated interchangeably in certain contexts.

setClassUnion(name = "AnyMatrixOrDataframe",
            members = c("matrix", "data.frame"))

#' Class \code{eQTLObject}
#'
#' The eQTLObject class is designed to store data related to eQTL analysis,
#' including data lists, result data frames, and additional metadata such as
#' classification, species, and grouping information.
#' @name eQTLObject-class
#' @slot rawData A gene expression dataframe, the row names represent gene IDs
#' and the column names represent cell IDs.
#' @slot filterData Gene expression matrix after normalizing.
#' @slot eQTLResult The result dataframe obtained the sc-eQTL results.
#' @slot biClassify Logical; whether to convert genotype encoding in snpMatrix
#' to 0, 1, and 2. \code{TRUE} indicates conversion; \code{FALSE} indicates no
#' conversion (default).
#' @slot species The species that the user wants to select, human or mouse.
#' @slot groupBy Options for cell grouping, users can choose celltype,
#' cellstatus, etc., depending on metadata.
#' @slot useModel model for fitting dataframe.
#'
#' @return An S4 object of class \code{eQTLObject}.
#' @export
setClass(Class = "eQTLObject",
        slots = c(rawData = "list",
                    filterData = "list",
                    eQTLResult = "AnyMatrixOrDataframe",
                    biClassify = "ANY",
                    species = "ANY",
                    groupBy = "data.frame",
                    useModel = "ANY"))


## Show Method for eQTLObject Class
##
## This method is to display information about an object of class eQTLObject.
## When called on an eQTLObject, it prints a descriptive message to the console
##
## @param object An S4 object of class eQTLObject.
## @importFrom methods setMethod
## @return information of eQTLObject
setMethod(f = "show", signature = "eQTLObject",
            definition = function(object) {
            cat("An object of class eQTLObject\n")
            cat("rawData:\n")
            cat(nrow(get_raw_data(object)[["rawExpMat"]]),
                "features across", ncol(get_raw_data(object)[["rawExpMat"]]),
                "samples in expressionMatrix\n")
            cat(nrow(get_raw_data(object)[["snpMat"]]),
                "SNPs across",
                ncol(get_raw_data(object)[["snpMat"]]),
                "samples in snpMatrix\n")
            if (length(get_filter_data(object)) != 0) {
                cat("filterData:\n")
                cat(nrow(get_filter_data(object)[["expMat"]]),
                    "features across",
                    ncol(get_filter_data(object)[["expMat"]]),
                    "samples in expressionMatrix\n")
                cat(nrow(get_filter_data(object)[["snpMat"]]),
                    "SNPs across",
                    ncol(get_filter_data(object)[["snpMat"]]),
                    "samples in snpMatrix")
            }
            if (length(get_model_info(object)) != 0) {
                cat("model used: ", get_model_info(object))
            }
    })


#' Create an eQTLObject for storing sc-eQTL analysis data.
#'
#' This function creates an S4 object to store genotype and expression data,
#' along with sample grouping and metadata for eQTL analysis.
#' @param snpMatrix A genotype matrix where each row is a snp and each column
#' is a cell. Encoding should be 0, 1, 2, 3.
#' @param genedata A gene expression matrix, or a Seurat object or
#' SingleCellExperiment object.
#' @param biClassify Logical; whether to convert genotype encoding in snpMatrix
#' to 0, 1, and 2. \code{TRUE} indicates conversion; \code{FALSE} indicates no
#' conversion (default).
#' @param species The species that the user wants to select, human or mouse.
#' @param group A column name in the metadata of the Seurat or
#' SingleCellExperiment object indicating cell groups (e.g., celltype or
#' condition).
#' @param ... other parameters
#' @importFrom SeuratObject GetAssayData
#' @importFrom methods is new
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assayNames assay
#' @return An S4 object of class \code{eQTLObject}.
#' @export
#'
#' @examples
#' data(SNPData)
#' data(GeneData)
#' eqtl <- createQTLObject(snpMatrix = SNPData,
#'                      genedata = GeneData,
#'                      biClassify = FALSE,
#'                      species = 'human',
#'                      group = NULL)
#'
createQTLObject <- function(snpMatrix,
                            genedata,
                            biClassify = FALSE,
                            species = NULL,
                            group = NULL,
                            ...) {
    expr <- extract_expr(genedata = genedata)
    raw_expressionMatrix <- expr[["raw_expressionMatrix"]]
    expressionMatrix <- expr[["expressionMatrix"]]

    check_raw_expr(expressionMatrix = raw_expressionMatrix)
    check_cell_names(expressionMatrix = raw_expressionMatrix,
                    snpMatrix = snpMatrix)

    expr_colnames <- colnames(raw_expressionMatrix)
    snp_colnames <- colnames(snpMatrix)

    # Extract cell grouping information from Seurat into metadata
    if (is.null(group)) {
        group_values <- rep("Matrix", length(expr_colnames))
        metadata <- data.frame(group = group_values)
        rownames(metadata) <- expr_colnames
    } else {
        if (is(genedata, "Seurat")) {
            metadata <- genedata[[group]]
            colnames(metadata) <- "group"
        } else if (is(genedata, "SingleCellExperiment")) {
            metadata <- colData(genedata, group)
            colnames(metadata) <- "group"
        } else {
            stop("Cannot find 'group'")
        }
    }
    raw_expressionMatrix <- as.matrix(raw_expressionMatrix)
    snpMatrix <- as.matrix(snpMatrix)

    assay <- new(Class = "eQTLObject",
                rawData = list(rawExpMat = raw_expressionMatrix,
                                normExpMat = expressionMatrix,
                                snpMat = snpMatrix),
                biClassify = biClassify,
                species = species,
                groupBy = metadata)
    return(assay)
}


# @rdname createQTLObject_internals
check_raw_expr <- function(expressionMatrix){
    if (sum(is.na(expressionMatrix)) > 0) {
        stop("NA detected in 'expressionMatrix'")
        gc()
    }
    if (sum(expressionMatrix < 0) > 0) {
        stop("Negative value detected in 'expressionMatrix'")
    }
    if (all(expressionMatrix == 0)) {
        stop("All elements of 'expressionMatrix' are zero")
    }
}


# @rdname createQTLObject_internals
check_cell_names <- function(expressionMatrix, snpMatrix){
    # Check if the cell name is duplicated
    if (anyDuplicated(colnames(expressionMatrix)) == 0 &
        anyDuplicated(colnames(snpMatrix)) == 0) {
        NULL
    } else {
        stop("There are duplicate values in the cell names.")
    }

    # Check whether the cell count and name of expressionMatrix match
    expr_colnames <- colnames(expressionMatrix)
    snp_colnames <- colnames(snpMatrix)
    if (all(expr_colnames %in% snp_colnames) &&
        all(snp_colnames %in% expr_colnames)) {
        NULL
    } else {
        stop("Column names do not match!")
    }
}


# @rdname createQTLObject_internals
extract_expr <- function(genedata){
    if (is(genedata, "Seurat")) {
        raw_expressionMatrix <- GetAssayData(genedata,
                                            assay = "RNA",
                                            layer = "counts")
        options(warn = -1)
        if (nrow(GetAssayData(genedata, assay = "RNA", layer = "data")) > 0) {
            expressionMatrix <- as.matrix(GetAssayData(genedata,
                                                        assay = "RNA",
                                                        layer = "data"))
        } else {
            message("No normalized counts in your object,
                you can achieve data normalization with `normalizeGene()`.")
            expressionMatrix <- NULL
        }
        options(warn = 1)
    } else if (is(genedata, "SingleCellExperiment")){
        if ("counts" %in% assayNames(genedata)) {
            raw_expressionMatrix <- assay(genedata, "counts")
        } else {
            stop("No raw count data found named 'counts'!
                Please check if it has been renamed.")
        }
        if ("normcounts" %in% assayNames(genedata)) {
            expressionMatrix <- assay(genedata, "normcounts")
        } else {
            message("No Normalized counts found named 'normcounts',
                you can achieve data normalization with `normalizeGene()`,
                or check if it is has been renamed.")
            expressionMatrix <- NULL
        }
    } else if (is.matrix(genedata) |
                is.data.frame(genedata) |
                is(genedata, "dgCMatrix")) {
        raw_expressionMatrix <- genedata
        expressionMatrix <- NULL
    } else {
        stop("Please enter the data in the correct format!")
    }
    return(list(raw_expressionMatrix = raw_expressionMatrix,
                expressionMatrix = expressionMatrix))
}
