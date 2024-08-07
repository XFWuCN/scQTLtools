#' visualizeQTL: Visualize the gene-snp pairs by group.
#' @param SNPid ID of SNP.
#' @param Geneid ID of Gene.
#' @param plottype Types of plot, one of  "QTLplot", "violin", "boxplot" or "histplot".
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param groupName Users can choose one or more than one single cell groups.
#' @param removeoutlier Whether identify and remove the outliers. Default by FALSE.
#'
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_violin geom_boxplot
#' @importFrom ggplot2 geom_point labs theme position_dodge scale_fill_brewer
#' @importFrom ggplot2 position_jitterdodge theme_bw element_text element_line
#' @importFrom ggplot2 unit element_blank geom_histogram facet_grid theme_minimal
#' @importFrom ggplot2 ggtitle guides guide_legend
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @return list
#' @export
#' @examples
#' data(testSNP)
#' data(testSNP)
#' eqtl <- createQTLObject(snpMatrix = testSNP,
#'                      genedata = testGene,
#'                      biClassify = FALSE,
#'                      species = 'human',
#'                      group = NULL)
#' eqtl <- normalizeGene(eqtl, method = "logNormalize")
#' eqtl <- filterGeneSNP(eqtl,
#'                       snp.number.of.cells.percent = 2,
#'                       expression.min = 0,
#'                       expression.number.of.cells.percent = 2
#'                       )
#' eqtl <- callQTL(eqtl,
#'                 gene_ids = NULL,
#'                 downstream = NULL,
#'                 upstream = NULL,
#'                 p.adjust.method = "bonferroni",
#'                 useModel = "poisson",
#'                 p.adjust.Threshold = 0.05,
#'                 logfc.threshold = 0.1
#'                 )
#' visualizeQTL(eqtl,
#'              SNPid = "1:632647",
#'              Geneid = "RPS27",
#'              groupName = NULL,
#'              plottype = "QTLplot",
#'              removeoutlier = FALSE
#'              )


visualizeQTL <- function(eQTLObject,
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
    split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == i]
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
        sample_gene = split_expressionMatrix[Geneid, ]
        non_zero_mask = sample_gene[sample_gene != 0]
        med = median(unlist(sample_no_zero))
        mad = mad(unlist(sample_no_zero))
        filter = sample_no_zero[colMeans(sample_no_zero) < med + 4 * mad]
        filter = as.data.frame(filter)
        ref_cells = intersect(ref_cells, rownames(filter))
        alt_cells = intersect(alt_cells, rownames(filter))
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
        sample_gene = split_expressionMatrix[Geneid,]
        sample_no_zero = sample_gene[sample_gene != 0]
        med = median(unlist(sample_no_zero))
        mad = mad(unlist(sample_no_zero))
        filter = sample_no_zero[sample_no_zero < med+4*mad]
        filter = as.data.frame(filter)
        AA_cells = intersect(AA_cells, rownames(filter))
        Aa_cells = intersect(Aa_cells, rownames(filter))
        aa_cells = intersect(aa_cells, rownames(filter))
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

  drawboxplot <- function(df, unique_group){
    ggplot(df, aes(x = factor(snp),
                   y = expression,
                   fill = factor(snp)))+
      geom_boxplot(alpha=0.3)+
      theme_bw() +
      scale_fill_brewer(palette="Dark2")+
      labs(title = unique_group, x = "", y = "Expression")+
      theme(axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.ticks = element_line(size=0.2, color="black"),
            axis.ticks.length = unit(0.2, "cm"),
            plot.title = element_text(hjust = 0.5, size = 14),
            legend.position = "none",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12))
  }

  drawviolinplot <- function(df, unique_group){
    ggplot(data = df,aes(x = snp,
                         y = expression ,
                         fill = factor(snp)))+
      scale_fill_manual(values = c("#D7AA36", "#D85356", "#94BBAD")) +
      geom_violin(alpha = 0.7, position = position_dodge(width = .75),
                  size = 0.8, color="black") +
      theme_bw() +
      labs(title = unique_group, y = "Expression", x = '') +
      theme(axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.ticks = element_line(size=0.2, color="black"),
            axis.ticks.length = unit(0.2, "cm"),
            plot.title = element_text(hjust = 0.5, size = 14),
            legend.position = "none",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12))
  }

  QTLplot <- function(df, unique_group){
    ggplot(data = df,aes(x = snp,
                         y = expression ,
                         fill = factor(snp)))+
      scale_fill_manual(values = c("#D7AA36", "#D85356", "#94BBAD")) +
      geom_violin(alpha = 0.7, position = position_dodge(width = .75),
                  size = 0.8, color="black") +
      geom_boxplot(notch = TRUE, outlier.size = -1,
                   color="black", lwd=0.5, alpha = 0.7) +
      geom_point(shape = 21, size=1.8, stroke=NA,
                 position = position_jitterdodge(jitter.width = 0.8),
                 alpha = 1.2) +
      theme_bw() +
      labs(title = unique_group, y = "Expression", x = "") +
      theme(axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.ticks = element_line(size=0.2, color="black"),
            axis.ticks.length = unit(0.2, "cm"),
            plot.title = element_text(hjust = 0.5, size = 14),
            legend.position = "none",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12))
  }

  drawhistplot <- function(df, unique_group){
    ggplot(df, aes(expression, fill = factor(snp)))+
      geom_histogram()+facet_grid(snp ~ ., margins=FALSE, scales="free_y") +
      scale_fill_brewer(palette = "Pastel1")+
      labs(title = unique_group,
           x = "Expression",
           y = "Count")+
      theme_minimal()+
      ggtitle(title)+
      theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
            axis.title.y = element_text(vjust = 0.5, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 14),
            legend.position = "none")+
      guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))
  }

  title_all <- paste("Plot", "of", Geneid , "and", SNPid )

  options(warn = -1)
  plot_list <- list()
  if(plottype=='QTLplot'){
    for(j in 1:length(unique_group)){
      group_split <- unique_group[j]
      title <- group_split
      df_split1 <- df_all$group == group_split
      df_split <- df_all[df_split1, ]
      plot <- QTLplot(df_split, unique_group[j])
      plot_list[[j]] <- plot
    }
  }else if(plottype=='violin'){
    for(j in 1:length(unique_group)){
      group_split <- unique_group[j]
      title <- group_split
      df_split1 <- df_all$group == group_split
      df_split <- df_all[df_split1, ]
      plot <- drawviolinplot(df_split, unique_group[j])
      plot_list[[j]] <- plot
    }
  }else if(plottype=='boxplot'){
    for(j in 1:length(unique_group)){
      group_split <- unique_group[j]
      title <- group_split
      df_split1 <- df_all$group == group_split
      df_split <- df_all[df_split1, ]
      plot <- drawboxplot(df_split, unique_group[j])
      plot_list[[j]] <- plot
    }
  }else if(plottype=='histplot'){
    for(j in 1:length(unique_group)){
      group_split <- unique_group[j]
      title <- group_split
      df_split1 <- df_all$group == group_split
      df_split <- df_all[df_split1, ]
      plot <- drawhistplot(df_split, unique_group[j])
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
  print(combined_plot)
  options(warn = 0)

}
