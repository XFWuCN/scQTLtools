#' Generate a boxplot of expression levels by SNP factor
#'
#' `draw_boxplot()` creates a boxplot visualizing expression levels across
#' different SNP factors in the dataframe. It uses ggplot2 to produce a plot
#' with customizable aesthetics for clarity and presentation.
#' @param df Data frames listed as gene expression data, genotype data,
#' and groups
#' @param unique_group name of unique group
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' set.seed(123)
#' counts_Ref <- rnorm(50, mean = 10, sd = 2)
#' counts_Alt <- rnorm(50, mean = 12, sd = 2)
#' i <- rep("GroupA", 100)
#' unique_group <- unique(i)
#' dataframe <- data.frame(
#' expression = c(counts_Ref, counts_Alt),
#' snp = c(rep("REF", length(counts_Ref)),
#'         rep("ALT", length(counts_Alt))), group = i)
#' dataframe$snp <- factor(dataframe$snp, levels = c("REF", "ALT"))
#' draw_boxplot(df = dataframe, unique_group = unique_group)
draw_boxplot <- function(df, unique_group) {
    snp <- df$snp
    ggplot(df, aes(x = factor(snp),
                    y = expression,
                    fill = factor(snp))) +
        geom_boxplot(alpha = 0.3) +
        theme_bw() +
        scale_fill_brewer(palette = "Dark2") +
        labs(title = unique_group, x = "", y = "Expression") +
        plots_theme_opts()
}
