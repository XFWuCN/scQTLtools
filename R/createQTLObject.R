# defines a new class union named 'AnyMatrixOrDataframe' This class
# union allows objects of class 'matrix', 'dgCMatrix', or
# 'data.frame' to be treated interchangeably in certain contexts.

setClassUnion(name = "AnyMatrixOrDataframe",
            members = c("matrix", "data.frame"))

#' Class 'eQTLObject'
#' The eQTLObject class is an R object designed to store data related to eQTL
#' analysis, encompassing data lists, result data frames, and layers for
#' biClassify, species, and grouping information.
#' @name eQTLObject-class
#' @slot rawData A gene expression dataframe, the row names represent gene IDs
#' and the column names represent cell IDs.
#' @slot filterData Gene expression matrix after normalizing.
#' @slot eQTLResult The result dataframe obtained the sc-eQTL results.
#' @slot biClassify The user chooses whether to convert the counting method of
#' the snpMatrix to 0, 1, 2, TRUE indicates conversion, and FALSE indicates no
#' conversion, default is no conversion.
#' @slot species The species that the user wants to select, human or mouse.
#' @slot groupBy Options for cell grouping, users can choose celltype,
#' cellstatus, etc., depending on metadata.
#' @slot useModel model for fitting dataframe.
#'
#' @return eQTLObject
#' @export
setClass(Class = "eQTLObject",
        slots = c(rawData = "list",
                    filterData = "list",
                    eQTLResult = "AnyMatrixOrDataframe",
                    biClassify = "ANY",
                    species = "ANY",
                    groupBy = "data.frame",
                    useModel = "ANY"))


#' Show Method for eQTLObject Class
#'
#' This method is to display information about an object of class eQTLObject.
#' When called on an eQTLObject, it prints a descriptive message to the console
#'
#' @param object An S4 object of class eQTLObject.
#' @importFrom methods setMethod
#' @return information of eQTLObject
#' @examples
#' data(testEQTL)
#' testEQTL
setMethod(f = "show", signature = "eQTLObject",
            definition = function(object) {
            cat("An object of class eQTLObject\n")
            cat("rawData:\n")
            cat(nrow(object@rawData$rawExpMat),
                "features across", ncol(object@rawData$rawExpMat),
                "samples in expressionMatrix\n")
            cat(nrow(object@rawData$snpMat), "SNPs across",
                ncol(object@rawData$snpMat), "samples in snpMatrix\n")
            if (length(object@filterData) != 0) {
                cat("filterData:\n")
                cat(nrow(object@filterData$expMat),
                    "features across", ncol(object@filterData$expMat),
                    "samples in expressionMatrix\n")
                cat(nrow(object@filterData$snpMat),
                    "SNPs across",
                    ncol(object@filterData$snpMat),
                    "samples in snpMatrix")
            }
            if (length(object@useModel) != 0) {
                cat("model used: ", object@useModel)
            }
            })


#' createObject:Create the eQTLObject.
#' We next create a S4 object. The object serves as a container that contains
#' both data (like the count matrix) and meta.data.
#'
#' @param snpMatrix A genotype matrix where each row is one variant and
#' each column is one sample, and the scoring method is 0/1/2/3.
#' @param genedata A gene expression matrix or a Seurat object, or a
#' SingleCellExperiment object.
#' @param biClassify The user chooses whether to convert the counting method of
#' the snpMatrix to 0/1/2, TRUE indicates conversion,
#' and FALSE indicates no conversion, default is no conversion.
#' @param species The species that the user wants to select, human or mouse.
#' @param group Provided by Seurat's meta.data, such as celltypes, cellstatus
#' and so on. By default, it is NULL.
#' @param ... other parameters
#' @importFrom SeuratObject GetAssayData
#' @importFrom methods is new
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assayNames
#' @return eQTLObject
#' @export
#'
#' @examples
#' data(testSNP)
#' data(testGene)
#' eqtl <- createQTLObject(snpMatrix = testSNP,
#'                      genedata = testGene,
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

    if (sum(is.na(raw_expressionMatrix)) > 0) {
        stop("NA detected in 'expressionMatrix'")
            gc()
    }
    if (sum(raw_expressionMatrix < 0) > 0) {
        stop("Negative value detected in 'expressionMatrix'")
    }
    if (all(raw_expressionMatrix == 0)) {
        stop("All elements of 'expressionMatrix' are zero")
    }

    # Check if the cell name is duplicated
    if (anyDuplicated(colnames(raw_expressionMatrix)) == 0 &
        anyDuplicated(colnames(snpMatrix)) == 0) {
        NULL
    } else {
        stop("There are duplicate values in the cell names.")
    }

    # Check whether the cell count and name of expressionMatrix match
    expr_colnames <- colnames(raw_expressionMatrix)
    snp_colnames <- colnames(snpMatrix)
    if (all(expr_colnames %in% snp_colnames) &&
        all(snp_colnames %in% expr_colnames)) {
        NULL
    } else {
        stop("Column names do not match!")
    }

    # Extract cell grouping information from Seurat into metadata
    if (is.null(group)) {
        group_values <- rep("Matrix", length(expr_colnames))
        metadata <- data.frame(group = group_values)
        rownames(metadata) <- expr_colnames
    } else {
        if (is(genedata, "Seurat")) {
            metadata <- genedata@meta.data[group]
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
