#' Visualize eQTL results by SNP–gene pair across groups
#'
#' @param SNPid ID of SNP.
#' @param Geneid ID of Gene.
#' @param plottype A \code{character} string specifying the type of plot to
#' generate. Must be one of \code{"QTLplot"}, \code{"violin"},
#' \code{"boxplot"}, or \code{"histplot"}.
#' @param eQTLObject An S4 object of class \code{eQTLObject}.
#' @param groupName A \code{character} vector specifying one or more cell group
#' names to include in the plot.
#' @param removeoutlier Logical; whether to identify and remove outliers.
#' Default is \code{FALSE}.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom stats median
#' @importFrom graphics title
#'
#' @return list
#' @export
#' @examples
#' data(EQTL_obj)
#' ## We have to call the eQTLs firstly using \code{\link{callQTL()}}.
#' eqtl <- callQTL(eQTLObject = EQTL_obj, useModel = "linear")
#' visualizeQTL(eQTLObject = eqtl,
#' SNPid = "1:632647",
#' Geneid = "RPS27",
#' groupName = NULL,
#' plottype = "QTLplot",
#' removeoutlier = FALSE)
visualizeQTL <- function(eQTLObject,
                        SNPid,
                        Geneid,
                        groupName = NULL,
                        plottype = 'QTLplot',
                        removeoutlier = FALSE) {
    eQTLresult <- get_result_info(eQTLObject)
    expressionMatrix <- get_filter_data(eQTLObject)[["expMat"]]
    snpMatrix <- get_filter_data(eQTLObject)[["snpMat"]]
    biClassify <- load_biclassify_info(eQTLObject)
    unique_group <- if (is.null(groupName)) {
        unique(load_group_info(eQTLObject)[["group"]])
    } else {
        groupName
    }
    process_group <- function(i) {
        split_cells <- rownames(load_group_info(eQTLObject)
        )[load_group_info(eQTLObject)[["group"]] == i]
        split_expressionMatrix <- expressionMatrix[, split_cells, drop = FALSE]
        snpMatrix_split <- snpMatrix[, split_cells, drop = FALSE]

        result <- eQTLresult[(eQTLresult$group == i) &
                            (eQTLresult$SNPid == SNPid) &
                            (eQTLresult$Geneid == Geneid), ]
        df <- process_classify(split_expressionMatrix, snpMatrix_split,
                                SNPid, Geneid, i, removeoutlier, biClassify)
        return(df)
    }
    df_all <- do.call(rbind, lapply(unique_group, process_group))
    title_all <- paste("Plot", "of", Geneid, "and", SNPid)
    options(warn = -1)
    plot_list <- lapply(unique_group, function(group_split) {
    df_split <- df_all[df_all$group == group_split, , drop = FALSE]
    switch(plottype,
            'QTLplot' = draw_QTLplot(df_split, group_split),
            'violin' = draw_violinplot(df_split, group_split),
            'boxplot' = draw_boxplot(df_split, group_split),
            'histplot' = draw_histplot(df_split, group_split),
            stop("Invalid plottype,
        Please choose from 'QTLplot', 'violin' , 'boxplot' or 'histplot'."))
    })
    combined_plot <- wrap_plots(plot_list)
    title_annotation <- plot_annotation(title = title_all,
                                        theme = theme(
                                        plot.title = element_text(
                                        size = 16, hjust = 0.5, vjust = 1)))
    combined_plot <- combined_plot + title_annotation
    options(warn = 0)
    return(combined_plot)
}


# @rdname visualizeQTL_internals
process_classify <- function(split_expressionMatrix,
                            snpMatrix_split,
                            SNPid,
                            Geneid,
                            group,
                            removeoutlier,
                            biClassify) {
    cell_groups <- get_cell_groups(snpMatrix_split, SNPid, biClassify)
    snp_labels <- if (biClassify) {
        c("REF", "ALT")
    } else {
        c("ref/ref", "ref/alt", "alt/alt")
    }
    if (removeoutlier) {
        do.call(remove_outliers,
            c(list(split_expressionMatrix, Geneid),
                cell_groups))
    }
    counts <- lapply(cell_groups,
                    get_counts,
                    expressionMatrix = split_expressionMatrix,
                    Geneid = Geneid)
    df <- data.frame(
        expression = unlist(counts),
        snp = rep(snp_labels, times = vapply(counts, length, integer(1))),
        group = group
    )
    df$snp <- factor(df$snp, levels = snp_labels)
    return(df)
}
