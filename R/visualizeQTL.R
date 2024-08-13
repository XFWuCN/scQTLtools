#' visualizeQTL: Visualize the gene-snp pairs by group.
#' @param SNPid ID of SNP.
#' @param Geneid ID of Gene.
#' @param plottype Types of plot,one of "QTLplot","violin","boxplot" or
#' "histplot".
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param groupName Users can choose one or more than one single cell groups.
#' @param removeoutlier Whether identify and remove the outliers.
#' Default by FALSE.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom stats median
#' @importFrom graphics title
#'
#' @return list
#' @export
#' @examples
#' data(testEQTL)
#' ## We have to call the eQTLs firstly using `callQTL()`.
#' eqtl <- callQTL(eQTLObject = testEQTL, useModel = "linear")
#' visualizeQTL(eQTLObject = eqtl,
#' SNPid = "1:632647",
#' Geneid = "RPS27",
#' groupName = NULL,
#' plottype = "QTLplot",
#' removeoutlier = FALSE)
visualizeQTL <- function(
    eQTLObject,
    SNPid,
    Geneid,
    groupName = NULL,
    plottype = 'QTLplot',
    removeoutlier = FALSE){

    eQTLresult <- eQTLObject@eQTLResult
    expressionMatrix <- eQTLObject@filterData$expMat
    snpMatrix <- eQTLObject@filterData$snpMat
    biClassify <- eQTLObject@biClassify

    if(is.null(groupName)){
    unique_group <- unique(eQTLObject@groupBy$group)
    }else{
        unique_group <- groupName
    }

    df_all <- data.frame()

    for(i in unique_group){
        split_cells <- rownames(eQTLObject@groupBy)[
        eQTLObject@groupBy$group == i]
        split_expressionMatrix <- expressionMatrix[, split_cells]
        snpMatrix_split <- snpMatrix[, split_cells]

    result <- eQTLresult[(eQTLresult$group == i)&
                        (eQTLresult$SNPid == SNPid)&
                        (eQTLresult$Geneid == Geneid),]

    if(biClassify == TRUE){

        result_split1 <- result$group == i
        result_split <- result[result_split1, ]

        adjust.pvalue <- result_split$adjusted_pvalue

        snpMatrix_split[snpMatrix_split == 3] <- 2

        ref_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid,] == 1]
        alt_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid,] == 2]

        if (removeoutlier) {
            remove_outliers(split_expressionMatrix,
                            Geneid,
                            ref_cells,
                            alt_cells)
        }

        counts_Ref <- unlist(split_expressionMatrix[Geneid, ref_cells])
        counts_Alt <- unlist(split_expressionMatrix[Geneid, alt_cells])

        # building data frame
        df <- data.frame(expression = c(counts_Ref,counts_Alt),
                        snp = c(rep("REF", length(counts_Ref)),
                                rep("ALT", length(counts_Alt))),
                        group = i)
        df$snp <- factor(df$snp, levels = c("REF", "ALT"))
        title <- paste(i)
        df_all <- rbind(df_all, df)

    }else if(biClassify == FALSE){

        result_split1 <- result$group == i
        result_split <- result[result_split1, ]

        adjust.pvalue <- result_split$adjusted_pvalue

        AA_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid,] == 1]
        Aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid,] == 3]
        aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid,] == 2]

        if (removeoutlier){
            remove_outliers(split_expressionMatrix,
                            Geneid,
                            AA_cells,
                            Aa_cells,
                            aa_cells)
        }

        counts_AA <- unlist(split_expressionMatrix[Geneid, AA_cells])
        counts_Aa <- unlist(split_expressionMatrix[Geneid, Aa_cells])
        counts_aa <- unlist(split_expressionMatrix[Geneid, aa_cells])

        # building data frame
        df <- data.frame(expression = c(counts_AA, counts_Aa, counts_aa),
                        snp = c(rep("ref/ref", length(counts_AA)),
                                rep("ref/alt", length(counts_Aa)),
                                rep("alt/alt", length(counts_aa))),
                        group = i)
        df$snp <- factor(df$snp, levels = c("ref/ref", "ref/alt", "alt/alt"))
        title <- paste(i)
        df_all <- rbind(df_all, df)

    }else{
        stop("biClassify can only be selected as 'TRUE' or 'FALSE'")
        }
    }

    title_all <- paste("Plot", "of", Geneid , "and", SNPid )

    options(warn = -1)
    plot_list <- list()
    if(plottype=='QTLplot'){
        for(j in seq_len(length(unique_group))){
            group_split <- unique_group[j]
            title <- group_split
            df_split <- df_all[df_all$group == group_split, ]
            snp <- df_split$snp
            plot <- draw_QTLplot(df_split, unique_group[j])
            plot_list[[j]] <- plot
        }
    }else if(plottype=='violin'){
        for(j in seq_len(length(unique_group))){
            group_split <- unique_group[j]
            title <- group_split
            df_split <- df_all[df_all$group == group_split, ]
            plot <- draw_violinplot(df_split, unique_group[j])
            plot_list[[j]] <- plot
        }
    }else if(plottype=='boxplot'){
    for(j in seq_len(length(unique_group))){
        group_split <- unique_group[j]
        title <- group_split
        df_split <- df_all[df_all$group == group_split, ]
        plot <- draw_boxplot(df_split, unique_group[j])
        plot_list[[j]] <- plot
    }
    }else if(plottype=='histplot'){
    for(j in seq_len(length(unique_group))){
        group_split <- unique_group[j]
        title <- group_split
        df_split <- df_all[df_all$group == group_split, ]
        plot <- draw_histplot(df_split, unique_group[j])
        plot_list[[j]] <- plot
    }
    }else{
    stop("Invalid plottype,
        Please choose from 'QTLplot', 'violin' , 'boxplot' or 'histplot'.")
    }
    combined_plot <- wrap_plots(plot_list)
    title_annotation <- plot_annotation(title = title_all,
                                        theme = theme(
                                        plot.title = element_text(
                                        size = 16, hjust = 0.5, vjust = 1)))
    combined_plot <- combined_plot + title_annotation
    return(combined_plot)
    options(warn = 0)
}
