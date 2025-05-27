#' Filter a data frame of SNP–gene pairs by absolute beta.
#'
#' @param result  A data frame that containing SNP–gene pairs and their beta
#' values (b-value).
#' @param logfcThreshold A numeric value specifying the minimum absolute
#' b-value to retain. Default is 0.1.
#'
#' @return A data frame containing only rows with absolute b-values above the
#' specified threshold.
#' @export
#'
#' @examples
#' example_result <- data.frame(
#'   gene = c("Gene1", "Gene2", "Gene3", "Gene4"),
#'   SNP = c("SNP1", "SNP2", "SNP3", "SNP4"),
#'   b = c(-2.5, 1.0, -0.5, 3.0))
#' logfcThreshold <- 0.1
#' filtered_result <- filter_by_abs_b(example_result, logfcThreshold)
filter_by_abs_b <- function(result, logfcThreshold) {
    result <- result %>%
        mutate(abs_b = abs(result[, "b"]))

    result <- result[result$abs_b >= logfcThreshold, , drop = FALSE]
    return(result)
}
