library(DESeq2)
library(limma)

# Load test data and create a mock eQTLObject
data(GeneData)
data(SNPData)
eqtl <- createQTLObject(snpMatrix = SNPData,
                        genedata = GeneData,
                        biClassify = FALSE,
                        species = 'human',
                        group = NULL)

# Test cases for normalizeGene function
test_that("normalizeGene function behaves as expected", {

  # Test logNormalize method
  eqtl_normalized <- normalizeGene(eqtl, method = "logNormalize")
  expect_true(is.null(eqtl_normalized@rawData$normExpMat) == FALSE)

  # Test CPM method
  eqtl_normalized <- normalizeGene(eqtl, method = "CPM")
  expect_true(is.null(eqtl_normalized@rawData$normExpMat) == FALSE)

  # Test TPM method
  eqtl_normalized <- normalizeGene(eqtl, method = "TPM")
  expect_true(is.null(eqtl_normalized@rawData$normExpMat) == FALSE)

  # Test DESeq method
  eqtl_normalized <- normalizeGene(eqtl, method = "DESeq")
  expect_true(is.null(eqtl_normalized@rawData$normExpMat) == FALSE)

  # Test limma method
  eqtl_normalized <- normalizeGene(eqtl, method = "limma")
  expect_true(is.null(eqtl_normalized@rawData$normExpMat) == FALSE)

  # Test invalid method
  expect_error(normalizeGene(eqtl, method = "invalidMethod"),
               "Invalid method.")

  # Test if input eQTLObject is normalized.
  eqtl_normalized <- normalizeGene(eqtl, method = "logNormalize")
  expect_false(identical(eqtl, eqtl_normalized))
})
