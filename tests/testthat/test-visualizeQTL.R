data(GeneData)
data(SNPData)
eqtl <- createQTLObject(snpMatrix = SNPData,
                        genedata = GeneData,
                        biClassify = FALSE,
                        species = 'human',
                        group = NULL)

eqtl <- normalizeGene(eqtl, method = "logNormalize")

eqtl <- filterGeneSNP(eQTLObject = eqtl,
                      snpNumOfCellsPercent = 2,
                      expressionMin = 0,
                      expressionNumOfCellsPercent = 2)

eqtl <- callQTL(eQTLObject = eqtl,
                gene_ids = NULL,
                downstream = NULL,
                upstream = NULL,
                pAdjustMethod = "bonferroni",
                useModel = "poisson",
                pAdjustThreshold = 1,
                logfcThreshold = 0)


test_that("visualizeQTL function behaves as expected", {

  # test QTLplot result
  plot1 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "QTLplot")
  expect_s3_class(plot1, "ggplot")

  # test violin result
  plot2 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "violin")
  expect_s3_class(plot2, "ggplot")

  # test boxplot result
  plot3 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "boxplot")
  expect_s3_class(plot3, "ggplot")

  # test histplot result
  plot4 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "histplot")
  expect_s3_class(plot4, "ggplot")

  # when specific group
  data(Seurat_obj)
  data(SNPData2)
  eqtl <- createQTLObject(snpMatrix = SNPData2,
                          genedata = Seurat_obj,
                          biClassify = FALSE,
                          species = 'human',
                          group = "celltype")
  eqtl <- normalizeGene(eqtl, method = "logNormalize")
  eqtl <- filterGeneSNP(eQTLObject = eqtl,
                        snpNumOfCellsPercent = 2,
                        expressionMin = 0,
                        expressionNumOfCellsPercent = 2)
  eqtl <- callQTL(eQTLObject = eqtl,
                  useModel = "linear",
                  pAdjustThreshold = 1,
                  logfcThreshold = 0)
  plot5 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", groupName = "GMP", plottype = "QTLplot")
  expect_s3_class(plot5, "ggplot")

  # test invalid plot types
  expect_error(visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "invalid_type"),
               "Invalid plottype,
        Please choose from 'QTLplot', 'violin' , 'boxplot' or 'histplot'.")
})
