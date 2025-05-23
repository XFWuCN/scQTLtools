#' Customized ggplot2 Theme for Plots
#'
#' @return A \code{ggplot2} theme object for styling plots.
#' @import ggplot2
#' @export
#'
#' @examples
#' library(ggplot2)
#' data <- data.frame(
#'   x = c("A", "B", "C", "D", "E"),
#'   y = c(10, 20, 30, 40, 50))
#' ggplot(data, aes(x, y)) +
#'   geom_point() +
#'   plots_theme_opts()
plots_theme_opts <- function(){
    theme(
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)
    )
}
