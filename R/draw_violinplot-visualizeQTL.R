#' Generate a violin plot of expression levels by SNP factor
#'
#' `draw_violinplot()` creates a violin plot visualizing expression levels
#' across different SNP factors in the dataframe. It uses ggplot2 to produce a
#' plot with customizable aesthetics for clarity and presentation.
#' @param df Data frames listed as gene expression data, genotype data,
#' and groups
#' @param unique_group name of unique group
#' @import ggplot2
#' @return ggplot
#' @export
#' @examples
#' set.seed(123)
#' counts_Ref <- rnorm(50, mean = 10, sd = 2)
#' counts_Alt <- rnorm(50, mean = 12, sd = 2)
#' i <- rep("GroupA", 100);unique_group <- unique(i)
#' dataframe <- data.frame(expression = c(counts_Ref, counts_Alt),
#'                         snp = c(rep("REF", length(counts_Ref)),
#'                         rep("ALT", length(counts_Alt))), group = i)
#' dataframe$snp <- factor(dataframe$snp, levels = c("REF", "ALT"))
#' draw_violinplot(df = dataframe, unique_group = unique_group)
draw_violinplot <- function(df, unique_group) {
    snp <- df$snp
    ggplot(data = df,
            aes(x = snp,
                y = expression,
                fill = factor(snp))) +
        scale_fill_manual(values = c("#D7AA36", "#D85356", "#94BBAD")) +
        geom_violin(alpha = 0.7, position = position_dodge(width = .75),
                    size = 0.8, color = "black") +
        theme_bw() +
        labs(title = unique_group, y = "Expression", x = "") +
        plots_theme_opts()
}
