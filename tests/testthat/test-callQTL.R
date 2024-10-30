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

test_that("Function returns correct output with valid gene_ids inputs", {
  result <- callQTL(eqtl, gene_ids = c("CNN2", "RNF113A"))
  unique_geneids <- unique(result@eQTLResult$Geneid)
  expect_length(unique_geneids, 2)
  expect_equal(unique_geneids, c("CNN2", "RNF113A"))
  expect_equal(rownames(result@filterData$snpMat),
               c("1:632647", "1:1203822",
                 "1:1407232", "1:1541864",
                 "1:1785578", "1:19623568",
                 "1:26281235", "1:26475513",
                 "1:28236165", "1:33013276",
                 "1:34865634", "1:43978026"))
})


test_that("Function handles invalid gene_ids gracefully", {
  expect_error(callQTL(eqtl, gene_ids = c("invalid_gene")),
    "The input gene_ids contain non-existent gene IDs.Please re-enter.")
})


test_that("Function handles upstream/downstream parameters correctly", {
  result <- callQTL(eqtl, downstream = -85000000, upstream = 20000000)
  expect_true("CCDC18" %in% result@eQTLResult$Geneid)
  boolean <- c("1:19623568", "1:26281235", "1:28236165", "1:43978026",
                "1:26475513", "1:33013276") %in% result@eQTLResult$SNPid
  expect_true(unique(boolean))
})


test_that("Function handles empty or invalid parameter", {
  eqtl@species <- NULL  # Simulate empty species
  expect_error(callQTL(eqtl, downstream = -85000000, upstream = 20000000),
               "The 'species' variable is NULL or empty.")

  eqtl@species <- "rat"  # Simulate invalid species
  expect_error(callQTL(eqtl, downstream = -85000000, upstream = 20000000),
               "Please enter 'human', 'mouse', 'worm' or 'phytozome'")

  eqtl@species <- "human"
  # invalid downstream
  expect_error(callQTL(eqtl, downstream = 85000000, upstream = 20000000),
               "downstream should be negative.")
  # do not enter upstream and downstream simultaneously
  expect_error(callQTL(eqtl, upstream = 20000000),
               "Please enter upstream and downstream simultaneously.")

  eqtl@filterData <- list()  # Simulate empty filterData
  expect_error(callQTL(eqtl), "Please filter the data first.")
})


test_that("Function returns correct output with valid SNP matrix rownames", {
  eqtl@filterData$snpMat <- eqtl@filterData$snpMat[1:5,]  # Simulate SNP matrix
  rownames(eqtl@filterData$snpMat) <- c("rs546",
                                        "rs549",
                                        "rs568",
                                        "rs665",
                                        "rs672")
  result <- callQTL(eqtl, downstream = -1e+6, upstream = 5e+6)
  expect_equal(result@eQTLResult$Geneid, "CCDC18")
  expect_equal(result@eQTLResult$SNPid, "rs546")
})

test_that("callQTL function behaves as expected when specific fitting model",{

  # test ZINB model
  result <- callQTL(
    eqtl,
    gene_ids = c("CNN2", "RNF113A", "TIGD2", "VWCE", "PLAU", "RPS27"),
    useModel = "zinb"
    )
  expect_true(is.null(result@eQTLResult) == FALSE)

  # test poisson model
  result <- callQTL(eqtl, useModel = "poisson")
  expect_true(is.null(result@eQTLResult) == FALSE)

  # test linear model
  result <- callQTL(eqtl, useModel = "linear")
  expect_true(is.null(result@eQTLResult) == FALSE)

  # test invalid method
  expect_error(callQTL(eqtl, useModel = "invalid method"),
    "Invalid model Please choose from 'zinb','poisson',or 'linear'.")
})


test_that("callQTL function behaves as expected when specific p value adjust
          method",{

  # test bonferroni method
  result <- callQTL(eqtl, pAdjustMethod = "bonferroni", useModel = "linear")
  expect_true(is.null(result@eQTLResult$adjusted_pvalue) == FALSE)

  # test holm method
  result <- callQTL(eqtl, pAdjustMethod = "holm", useModel = "linear")
  expect_true(is.null(result@eQTLResult$adjusted_pvalue) == FALSE)

  # test hochberg method
  result <- callQTL(eqtl, pAdjustMethod = "hochberg", useModel = "linear")
  expect_true(is.null(result@eQTLResult$adjusted_pvalue) == FALSE)

  # test hommel method
  result <- callQTL(eqtl, pAdjustMethod = "hommel", useModel = "linear")
  expect_true(is.null(result@eQTLResult$adjusted_pvalue) == FALSE)

  # test BH method
  result <- callQTL(eqtl, pAdjustMethod = "BH", useModel = "linear")
  expect_true(is.null(result@eQTLResult$adjusted_pvalue) == FALSE)

})
