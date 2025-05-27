library(methods)

# Load test data and create a mock eQTLObject
data(GeneData)
data(SNPData)
data(Seurat_obj)
data(SNPData2)

test_that("createQTLObject handles matrix input correctly", {
  eqtl <- createQTLObject(snpMatrix = SNPData, genedata = GeneData)
  expect_true(inherits(eqtl, "eQTLObject"))
  expect_equal(nrow(eqtl@rawData$rawExpMat), 100)
  expect_equal(ncol(eqtl@rawData$rawExpMat), 2705)
  expect_equal(nrow(eqtl@rawData$snpMat), 1000)
  expect_equal(ncol(eqtl@rawData$snpMat), 2705)
})

test_that("createQTLObject handles Seurat object input correctly", {
  eqtl <- createQTLObject(snpMatrix = SNPData2, genedata = Seurat_obj)
  expect_true(inherits(eqtl, "eQTLObject"))
  expect_equal(nrow(eqtl@rawData$rawExpMat), 100)
  expect_equal(ncol(eqtl@rawData$rawExpMat), 500)
  expect_equal(nrow(eqtl@rawData$snpMat), 500)
  expect_equal(ncol(eqtl@rawData$snpMat), 500)
  expect_equal(nrow(eqtl@groupBy), 500)
})

test_that("createQTLObject handles SingelCellExperiment object input correctly", {
  library(SingleCellExperiment)
  sce <- SingleCellExperiment(assays = list(counts = GeneData))
  eqtl <- createQTLObject(snpMatrix = SNPData, genedata = sce)
  expect_equal(nrow(eqtl@rawData$rawExpMat), 100)
  expect_equal(ncol(eqtl@rawData$rawExpMat), 2705)
  expect_equal(nrow(eqtl@rawData$snpMat), 1000)
  expect_equal(ncol(eqtl@rawData$snpMat), 2705)
})

test_that("createQTLObject handles incorrect input gracefully", {
  expect_error(createQTLObject(snpMatrix = SNPData, genedata = "invalid_input"))
})
