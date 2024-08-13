data(testGene)
data(testSNP)
eqtl <- createQTLObject(snpMatrix = testSNP,
                        genedata = testGene,
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
                pAdjustThreshold = 0.05,
                logfcThreshold = 0.1)


test_that("visualizeQTL function behaves as expected", {

  # test QTLplot result
  plot1 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "QTLplot")
  expect_true(is.list(plot1))

  # test violin result
  plot2 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "violin")
  expect_true(is.list(plot2))

  # test boxplot result
  plot3 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "boxplot")
  expect_true(is.list(plot3))

  # test histplot result
  plot4 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "histplot")
  expect_true(is.list(plot4))

  # when specific group
  data(testSeurat)
  data("testSNP2")
  eqtl <- createQTLObject(snpMatrix = testSNP2,
                          genedata = testSeurat,
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
                  pAdjustThreshold = 0.05,
                  logfcThreshold = 0.025)
  plot5 <- visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", groupName = "GMP", plottype = "QTLplot")
  expect_true(is.list(plot5))

  # test invalid plot types
  expect_error(visualizeQTL(eqtl, SNPid = "1:632647", Geneid = "RPS27", plottype = "invalid_type"),
               "Invalid plottype,
        Please choose from 'QTLplot', 'violin' , 'boxplot' or 'histplot'.")
})
