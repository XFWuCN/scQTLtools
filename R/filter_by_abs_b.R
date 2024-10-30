#' Filters data frame by absolute b-values, returning rows meeting or
#' exceeding a threshold.
#'
#' @param result  Dataframe that contains gene-SNP pairs' information.
#' @param logfcThreshold  Represents the minimum beta threshold for fitting
#' SNP-Gene pairs. Default by 0.1.
#'
#' @return A dataframe filtered by absolute b-values.
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
