library(progress)

# Load test data and create a mock eQTLObject
data(testGene)
data(testSNP)
eqtl <- createQTLObject(snpMatrix = testSNP,
                        genedata = testGene,
                        biClassify = FALSE,
                        species = 'human',
                        group = NULL)

# Unit tests for filterGeneSNP function
test_that("Normalization must be performed before executing the filterGeneSNP function.",{
  # No valid expression data (should throw an error)
  expect_error(filterGeneSNP(eqtl), "Please normalize the raw expression data first.")
})


test_that("filterGeneSNP correctly filters SNP and expression matrices", {

  # Normalize the eQTLObject (mocking the normalization function)
  eqtl <- normalizeGene(eqtl, method = "logNormalize")

  # Test case 1: Default parameters
  eqtl_filtered <- filterGeneSNP(eqtl)
  expect_true(is.matrix(eqtl_filtered@filterData$expMat))
  expect_true(is.matrix(eqtl_filtered@filterData$snpMat))

  # Test case 2: Adjusted parameters
  eqtl_filtered <- filterGeneSNP(eqtl,
                                 snpNumOfCellsPercent = 3,
                                 expressionMin = 0.5,
                                 expressionNumOfCellsPercent = 3)
  expect_true(is.matrix(eqtl_filtered@filterData$expMat))
  expect_true(is.matrix(eqtl_filtered@filterData$snpMat))

})
