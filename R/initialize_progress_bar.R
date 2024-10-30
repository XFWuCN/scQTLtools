#' Progress Bar for Model Analysis.
#'
#' This function initializes a progress bar for use in the `linearModel`
#' , `poissonModel` and `zinbModel` function. It is designed to provide
#' feedback on the progress of the analysis by displaying the current step and
#' a percentage completion.
#'
#' @param total The total number of steps or iterations for which the progress
#' bar will be updated.
#' @param k A label or identifier for the specific group or iteration for which
#' the progress bar is being initialized.
#'
#' @return A `progress_bar` object from the `progress` package,
#' which is used to track and display the progress.
#'
#' @export
#'
#' @examples
#' unique_group <- c("CMP", "GMP")
#' total_snp_count <- 10  # assume each group have 100 SNP.
#' pb_model <- lapply(unique_group, function(k) {
#'     pb <- initialize_progress_bar(total = total_snp_count, k)
#'     for (i in seq_len(total_snp_count)) {
#'         Sys.sleep(0.1)  # assume progress time
#'         pb$tick()  # update pb
#'     }
#' })
initialize_progress_bar <- function(total, k) {
    message(k, ":")
    message("0%   10   20   30   40   50   60   70   80   90   100%")
    message("[----|----|----|----|----|----|----|----|----|----|")
    pb <- progress_bar$new(
        total = total,
        format = "[:bar]",
        clear = FALSE,
        width = 51
        )
    return(pb)
}
