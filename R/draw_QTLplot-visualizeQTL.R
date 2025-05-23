#' Create a combined violin-box-scatter plot.
#'
#' \code{draw_QTLplot()} creates a composite plot that overlays violin plots,
#' boxplots, and scatter points to illustrate the distribution and variability
#' of gene expression across SNP groups.
#' @param df A data frame containing gene expression values, SNP genotypes,
#' and group labels.
#' @param unique_group A character string indicating the unique group name.
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
#' draw_QTLplot(df = dataframe, unique_group = unique_group)
draw_QTLplot <- function(df, unique_group) {
    snp <- df$snp
    ggplot(data = df, aes(
        x = snp,
        y = expression,
        fill = factor(snp))) +
        scale_fill_manual(values = c("#D7AA36", "#D85356", "#94BBAD")) +
        geom_violin(alpha = 0.7, position = position_dodge(width = .75),
                    size = 0.8, color = "black") +
        geom_boxplot(notch = FALSE, outlier.size = -1, color = "black",
                    lwd = 0.5, alpha = 0.7) +
        geom_point(position = position_jitterdodge(jitter.width = 0.8),
                    stroke = NA, shape = 21, size = 1.8, alpha = 1.2) +
        theme_bw() +
        labs(title = unique_group, y = "Expression", x = "") +
        plots_theme_opts()
}
