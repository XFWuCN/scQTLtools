library(methods)

# Load test data and create a mock eQTLObject
data(testGene)
data(testSNP)
data(testSeurat)
data(testSNP2)

test_that("createQTLObject handles matrix input correctly", {
  eqtl <- createQTLObject(snpMatrix = testSNP, genedata = testGene)
  expect_true(inherits(eqtl, "eQTLObject"))
  expect_equal(nrow(eqtl@rawData$rawExpMat), 100)
  expect_equal(ncol(eqtl@rawData$rawExpMat), 2705)
  expect_equal(nrow(eqtl@rawData$snpMat), 1000)
  expect_equal(ncol(eqtl@rawData$snpMat), 2705)
})

test_that("createQTLObject handles Seurat object input correctly", {
  eqtl <- createQTLObject(snpMatrix = testSNP2, genedata = testSeurat)
  expect_true(inherits(eqtl, "eQTLObject"))
  expect_equal(nrow(eqtl@rawData$rawExpMat), 100)
  expect_equal(ncol(eqtl@rawData$rawExpMat), 500)
  expect_equal(nrow(eqtl@rawData$snpMat), 500)
  expect_equal(ncol(eqtl@rawData$snpMat), 500)
  expect_equal(nrow(eqtl@groupBy), 500)
})

test_that("createQTLObject handles incorrect input gracefully", {
  expect_error(createQTLObject(snpMatrix = testSNP, genedata = "invalid_input"))
})
