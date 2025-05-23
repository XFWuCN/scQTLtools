#' Generate a histogram of gene expression by SNP genotype.
#'
#' \code{draw_histplot()} shows the distribution of expression values for each
#' SNP genotype using histograms.
#' @param df A data frame containing gene expression values, SNP genotypes,
#' and group labels.
#' @param unique_group A character string indicating the unique group name.
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' set.seed(123)
#' counts_Ref <- rnorm(50, mean = 10, sd = 2)
#' counts_Alt <- rnorm(50, mean = 12, sd = 2)
#' i <- rep("GroupA", 100);unique_group <- unique(i)
#' dataframe <- data.frame(expression = c(counts_Ref, counts_Alt),
#'                         snp = c(rep("REF", length(counts_Ref)),
#'                         rep("ALT", length(counts_Alt))), group = i)
#' dataframe$snp <- factor(dataframe$snp, levels = c("REF", "ALT"))
#' draw_histplot(df = dataframe, unique_group = unique_group)
draw_histplot <- function(df, unique_group){
    df <- df; expression <- df$expression; snp <- df$snp
    ggplot(df, aes(expression, fill = snp))+
        geom_histogram()+
        facet_grid(snp ~ ., margins = FALSE, scales = "free_y") +
        scale_fill_brewer(palette = "Pastel1")+
        labs(title = unique_group,
        x = "Expression",
        y = "Count")+
    theme_minimal()+
    theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.title.y = element_text(vjust = 0.5, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "none")+
    guides(fill = guide_legend(title = NULL),
            color = guide_legend(title = NULL))
}
