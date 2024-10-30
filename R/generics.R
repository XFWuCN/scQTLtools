#' Generic to access eQTLObject raw data
#' @param x A eQTLObject object.
#'
#' @examples
#'   data(testEQTL)
#'   get_raw_data(testEQTL)
#'
#' @return raw data matrix.
#'
#' @export
setGeneric("get_raw_data", function(x) standardGeneric("get_raw_data"))

#' Method to access eQTLObject raw data
#' @param x A eQTLObject object.
#'
#' @return raw data matrix.
#'
#' @export
setMethod("get_raw_data", "eQTLObject", function(x) {
    value <- x@rawData
    return(value)
})

#' Generic to set eQTLObject raw data
#' @param x A eQTLObject object.
#' @param value The raw data.
#' @param name The matrix named 'name' is stored under the 'rawData' slot as an
#' element within its list.
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

#' Method to set eQTLObject raw data
#' @param x A eQTLObject object.
#' @param value The raw data.
#' @param name The matrix named 'name' is stored under the 'rawData' slot as an
#' element within its list.
#'
#' @examples
#'   data(testEQTL)
#'   data123 <- matrix(0, nrow = 3, ncol = 3)
#'   set_raw_data(testEQTL, data123, "rawExpMat")
#'
#' @return eQTLObject.
#'
#' @export
setMethod("set_raw_data", "eQTLObject", function(x, value, name) {
    x@rawData[[name]] <- value
    return(x)
})

#' Generic to access eQTLObject filter data
#' @param x A eQTLObject object.
#'
#' @examples
#'   data(testEQTL)
#'   get_filter_data(testEQTL)
#'
#' @return filtered matrices.
#'
#' @export
setGeneric("get_filter_data", function(x) standardGeneric("get_filter_data"))

#' Method to access eQTLObject filter data
#' @param x A eQTLObject object.
#'
#' @return filtered matrices.
#'
#' @export
setMethod("get_filter_data", "eQTLObject", function(x) {
    value <- x@filterData
    return(value)
})

#' Generic to set eQTLObject filter data
#' @param x A eQTLObject object.
#' @param value The filtered data.
#' @param name The matrix named 'name' is stored under the 'filterData' slot as
#' an element within its list.
#'
#' @examples
#'   data(testEQTL)
#'   data123 <- matrix(0, nrow = 3, ncol = 3)
#'   set_filter_data(testEQTL, data123, "expMat")
#'
#' @return eQTLObject.
#'
#' @export
setGeneric("set_filter_data", function(x, value, name)
    standardGeneric("set_filter_data"))

#' Method to set eQTLObject filter data
#' @param x A eQTLObject object.
#' @param value The filtered data.
#' @param name The matrix named 'name' is stored under the 'filterData' slot as
#' an element within its list.
#'
#' @examples
#'   data(testEQTL)
#'   data123 <- matrix(0, nrow = 3, ncol = 3)
#'   set_filter_data(testEQTL, data123, "expMat")
#'
#' @return eQTLObject.
#'
#' @export
setMethod("set_filter_data", "eQTLObject", function(x, value, name) {
    x@filterData[[name]] <- value
    return(x)
})

#' Generic to access the result of identifying eQTLs from scRNA-seq data
#' @param x A eQTLObject object.
#'
#' @examples
#'   data(testEQTL)
#'   get_result_info(testEQTL)
#'
#' @return A dataframe.
#'
#' @export
setGeneric("get_result_info", function(x) standardGeneric("get_result_info"))

#' Method to access the result of identifying eQTLs from scRNA-seq data
#' @param x A eQTLObject object.
#'
#' @return A dataframe.
#'
#' @export
setMethod("get_result_info", "eQTLObject", function(x) {
    value <- x@eQTLResult
    return(value)
})

#' Generic to set the result of identifying eQTLs from scRNA-seq data
#' @param x A eQTLObject object.
#' @param value A dataframe, each row describes eQTL discovering result of a
#' SNP-Gene pair.
#'
#' @examples
#'   data(testEQTL)
#'   result <- matrix(0, nrow = 3, ncol = 3)
#'   set_result_info(testEQTL, result)
#'
#' @return eQTLObject.
#'
#' @export
setGeneric("set_result_info", function(x, value)
    standardGeneric("set_result_info"))

#' Method to set the result of identifying eQTLs from scRNA-seq data
#' @param x A eQTLObject object.
#' @param value A dataframe, each row describes eQTL discovering result of a
#' SNP-Gene pair.
#'
#' @examples
#'   data(testEQTL)
#'   result <- matrix(0, nrow = 3, ncol = 3)
#'   set_result_info(testEQTL, result)
#'
#' @return eQTLObject.
#'
#' @export
setMethod("set_result_info", "eQTLObject", function(x, value) {
    x@eQTLResult <- value
    return(x)
})

#' Generic to access eQTLObject biclassify information
#' @param x A eQTLObject object.
#'
#' @examples
#'   data(testEQTL)
#'   load_biclassify_info(testEQTL)
#'
#' @return biclassify information of eQTLObject.
#'
#' @export
setGeneric("load_biclassify_info", function(x)
    standardGeneric("load_biclassify_info"))

#' Method to access eQTLObject biclassify information
#' @param x A eQTLObject object.
#'
#' @return biclassify information of eQTLObject.
#'
#' @export
setMethod("load_biclassify_info", "eQTLObject", function(x) {
    value <- x@biClassify
    return(value)
})

#' Generic to access eQTLObject species information
#' @param x A eQTLObject object.
#'
#' @examples
#'   data(testEQTL)
#'   load_species_info(testEQTL)
#'
#' @return species information of eQTLObject.
#'
#' @export
setGeneric("load_species_info", function(x)
    standardGeneric("load_species_info"))

#' Method to access eQTLObject species information
#' @param x A eQTLObject object.
#'
#' @return species information of eQTLObject.
#'
#' @export
setMethod("load_species_info", "eQTLObject", function(x) {
    value <- x@species
    return(value)
})

#' Generic to access eQTLObject cell grouping information
#' @param x A eQTLObject object.
#'
#' @examples
#'   data(testEQTL)
#'   load_group_info(testEQTL)
#'
#' @return A dataframe.
#'
#' @export
setGeneric("load_group_info", function(x) standardGeneric("load_group_info"))

#' Method to access eQTLObject cell grouping information
#' @param x A eQTLObject object.
#'
#' @return A dataframe.
#'
#' @export
setMethod("load_group_info", "eQTLObject", function(x) {
    value <- x@groupBy
    return(value)
})

#' Generic to access eQTLObject used model information
#' @param x A eQTLObject object.
#'
#' @examples
#'   data(testEQTL)
#'   get_model_info(testEQTL)
#'
#' @return used model information of eQTLObject.
#'
#' @export
setGeneric("get_model_info", function(x) standardGeneric("get_model_info"))

#' Method to access eQTLObject used model information
#' @param x A eQTLObject object.
#'
#' @return used model information of eQTLObject.
#'
#' @export
setMethod("get_model_info", "eQTLObject", function(x) {
    value <- x@useModel
    return(value)
})

#' Generic to set eQTLObject used model information
#' @param x A eQTLObject object.
#' @param value The used model information to set to eQTLObject.
#'
#' @examples
#'   data(testEQTL)
#'   useModel <- "zinb"
#'   set_model_info(testEQTL, useModel)
#'
#' @return eQTLObject.
#'
#' @export
setGeneric("set_model_info", function(x, value)
    standardGeneric("set_model_info"))

#' Method to set eQTLObject used model information
#' @param x A eQTLObject object.
#' @param value The used model information to set to eQTLObject.
#'
#' @examples
#'   data(testEQTL)
#'   useModel <- "zinb"
#'   set_model_info(testEQTL, useModel)
#'
#' @return eQTLObject.
#'
#' @export
setMethod("set_model_info", "eQTLObject", function(x, value) {
    x@useModel <- value
    return(x)
})
