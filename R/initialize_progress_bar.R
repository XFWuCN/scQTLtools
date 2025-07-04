#' Progress Bar for Model Analysis.
#'
#' Initializes a progress bar to track the computation progress during model
#' fitting procedures such as \code{linearModel}, \code{poissonModel}, and
#' \code{zinbModel}. This helps users monitor the status of long-running
#' analyses.
#'
#' @param total Integer. The total number of iterations or steps to be
#' completed in the analysis.
#' @param k Character. An identifier or label for the specific group or subset
#' being analyzed, used to annotate progress messages.
#'
#' @return An object of class \code{progress_bar} from the \pkg{progress}
#' package, which can be updated using the \code{$tick()} method.
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
