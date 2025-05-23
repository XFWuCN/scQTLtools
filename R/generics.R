#' Access raw data stored in an eQTLObject.
#' @param x An \code{eQTLObject}.
#'
#' @examples
#'   data(testEQTL)
#'   get_raw_data(testEQTL)
#'
#' @return raw data matrix.
#'
#' @export
setGeneric("get_raw_data", function(x) standardGeneric("get_raw_data"))

#' Method to access eQTLObject raw data.
#' @param x An \code{eQTLObject}.
#'
#' @return raw data matrix.
#'
#' @rdname get_raw_data
#' @export
setMethod("get_raw_data", "eQTLObject", function(x) {
    value <- x@rawData
    return(value)
})


#' Set raw data in an eQTLObject.
#' @param x An \code{eQTLObject}.
#' @param value A matrix to be stored as raw data.
#' @param name A character string indicating the key under which the matrix is
#' stored in the \code{rawData} list.
#'
#' @examples
#'   data(testEQTL)
#'   data123 <- matrix(0, nrow = 3, ncol = 3)
#'   set_raw_data(testEQTL, data123, "rawExpMat")
#'
#' @return eQTLObject.
#'
#' @export
setGeneric("set_raw_data", function(x, value, name)
    standardGeneric("set_raw_data"))

#' Method to set eQTLObject raw data.
#' @param x An \code{eQTLObject}.
#' @param value The raw data.
#' @param name A character string indicating the key under which the matrix is
#' stored in \code{rawData}.
#'
#' @examples
#'   data(testEQTL)
#'   data123 <- matrix(0, nrow = 3, ncol = 3)
#'   set_raw_data(testEQTL, data123, "rawExpMat")
#'
#' @return An updated \code{eQTLObject}.
#' @rdname set_raw_data
#' @export
setMethod("set_raw_data", "eQTLObject", function(x, value, name) {
    x@rawData[[name]] <- value
    return(x)
})


#' Access filtered data stored in an eQTLObject.
#' @param x An \code{eQTLObject}.
#'
#' @examples
#'   data(testEQTL)
#'   get_filter_data(testEQTL)
#'
#' @return filtered matrices.
#'
#' @export
setGeneric("get_filter_data", function(x) standardGeneric("get_filter_data"))

#' Method to access eQTLObject filter data.
#' @param x An \code{eQTLObject}.
#'
#' @return Filtered matrices.
#' @rdname get_filter_data
#' @export
setMethod("get_filter_data", "eQTLObject", function(x) {
    value <- x@filterData
    return(value)
})


#' Set filtered data in an eQTLObject.
#' @param x An \code{eQTLObject}.
#' @param value A matrix to be stored as filtered data.
#' @param name A character string indicating the key under which the matrix is
#' stored in the \code{filterData} list.
#'
#' @examples
#'   data(testEQTL)
#'   data123 <- matrix(0, nrow = 3, ncol = 3)
#'   set_filter_data(testEQTL, data123, "expMat")
#'
#' @return An updated \code{eQTLObject}.
#'
#' @export
setGeneric("set_filter_data", function(x, value, name)
    standardGeneric("set_filter_data"))

#' Method to set eQTLObject filter data.
#' @param x An \code{eQTLObject}.
#' @param value A matrix to be stored as filtered data.
#' @param name A character string indicating the key under which the matrix is
#' stored in \code{filterData}.
#'
#' @examples
#'   data(testEQTL)
#'   data123 <- matrix(0, nrow = 3, ncol = 3)
#'   set_filter_data(testEQTL, data123, "expMat")
#'
#' @return An updated \code{eQTLObject}.
#' @rdname set_filter_data
#' @export
setMethod("set_filter_data", "eQTLObject", function(x, value, name) {
    x@filterData[[name]] <- value
    return(x)
})


#' Access eQTLs results from an eQTLObject.
#' @param x An \code{eQTLObject}.
#'
#' @examples
#'   data(testEQTL)
#'   get_result_info(testEQTL)
#'
#' @return A data frame where each row corresponds to an identified gene-SNP
#' pair.
#'
#' @export
setGeneric("get_result_info", function(x) standardGeneric("get_result_info"))

#' Method to access the result of identifying eQTLs.
#' @param x An \code{eQTLObject}.
#'
#' @return A data frame with eQTL results.
#' @rdname get_result_info
#' @export
setMethod("get_result_info", "eQTLObject", function(x) {
    value <- x@eQTLResult
    return(value)
})


#' Set eQTL results in an eQTLObject.
#' @param x An \code{eQTLObject}.
#' @param value A data frame in which each row describes the result for a
#' gene-SNP pair.
#'
#' @examples
#'   data(testEQTL)
#'   result <- data.frame(0, nrow = 3, ncol = 3)
#'   set_result_info(testEQTL, result)
#'
#' @return An updated \code{eQTLObject}.
#'
#' @export
setGeneric("set_result_info", function(x, value)
    standardGeneric("set_result_info"))

#' Method to set the result of identifying eQTLs from scRNA-seq data.
#' @param x An \code{eQTLObject}.
#' @param value A data frame where each row corresponds to a gene-SNP pair.
#'
#' @examples
#'   data(testEQTL)
#'   result <- matrix(0, nrow = 3, ncol = 3)
#'   set_result_info(testEQTL, result)
#'
#' @return An updated \code{eQTLObject}.
#' @rdname set_result_info
#' @export
setMethod("set_result_info", "eQTLObject", function(x, value) {
    x@eQTLResult <- value
    return(x)
})


#' Access biclassification information from an eQTLObject.
#' @param x An \code{eQTLObject}.
#'
#' @examples
#'   data(testEQTL)
#'   load_biclassify_info(testEQTL)
#' @return A character or list containing biclassification information.
#'
#' @export
setGeneric("load_biclassify_info", function(x)
    standardGeneric("load_biclassify_info"))

#' Method to access eQTLObject biclassify information.
#' @param x An \code{eQTLObject}.
#'
#' @return biclassify information of eQTLObject.
#' @rdname load_biclassify_info
#' @export
setMethod("load_biclassify_info", "eQTLObject", function(x) {
    value <- x@biClassify
    return(value)
})


#' Access species information from an eQTLObject.
#'
#' @param x An \code{eQTLObject}.
#'
#' @examples
#'   data(testEQTL)
#'   load_species_info(testEQTL)
#'
#' @return A character string indicating the species.
#'
#' @export
setGeneric("load_species_info", function(x)
    standardGeneric("load_species_info"))

#' Method to access eQTLObject species information.
#' @param x An \code{eQTLObject}.
#'
#' @return species information of eQTLObject.
#' @rdname load_species_info
#' @export
setMethod("load_species_info", "eQTLObject", function(x) {
    value <- x@species
    return(value)
})


#' Access cell grouping information from an eQTLObject
#'
#' @param x An \code{eQTLObject}.
#'
#' @examples
#'   data(testEQTL)
#'   load_group_info(testEQTL)
#'
#' @return A data frame with grouping information.
#'
#' @export
setGeneric("load_group_info", function(x) standardGeneric("load_group_info"))

#' Method to access eQTLObject cell grouping information.
#' @param x An \code{eQTLObject}.
#'
#' @return A data frame with cell group assignments.
#' @rdname load_group_info
#' @export
setMethod("load_group_info", "eQTLObject", function(x) {
    value <- x@groupBy
    return(value)
})


#' Access model specification from an eQTLObject
#'
#' @param x An \code{eQTLObject}.
#'
#' @examples
#'   data(testEQTL)
#'   get_model_info(testEQTL)
#'
#' @return A character string indicating the model.
#'
#' @export
setGeneric("get_model_info", function(x) standardGeneric("get_model_info"))

#' Method to access eQTLObject used model information.
#' @param x An \code{eQTLObject}.
#'
#' @return used model information of eQTLObject.
#' @rdname get_model_info
#' @export
setMethod("get_model_info", "eQTLObject", function(x) {
    value <- x@useModel
    return(value)
})


#' Set model specification in an eQTLObject
#'
#' @param x An \code{eQTLObject}.
#' @param value A character string indicating the model (e.g., \code{"zinb"},
#' \code{"poisson"}, or \code{"linear"}).
#'
#' @examples
#'   data(testEQTL)
#'   useModel <- "zinb"
#'   set_model_info(testEQTL, useModel)
#'
#' @return An updated \code{eQTLObject} with the new model specification.
#'
#' @export
setGeneric("set_model_info", function(x, value)
    standardGeneric("set_model_info"))

#' Method to set eQTLObject used model information.
#' @param x An \code{eQTLObject}.
#' @param value A character string specifying the model (e.g., "zinb",
#' "poisson", "linear").
#'
#' @examples
#'   data(testEQTL)
#'   useModel <- "zinb"
#'   set_model_info(testEQTL, useModel)
#'
#' @return An updated \code{eQTLObject}.
#' @rdname set_model_info
#' @export
setMethod("set_model_info", "eQTLObject", function(x, value) {
    x@useModel <- value
    return(x)
})
