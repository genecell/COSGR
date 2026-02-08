context("Identification of Marker Genes by COSG")

library(COSG)
library(SeuratObject)

# Load test data from SeuratObject package (compatible with v5)
data("pbmc_small", package = "SeuratObject")

test_that("COSG works with default parameters", {
  expect_error(cosg(object = pbmc_small), NA)
})

test_that("COSG works with layer parameter (Seurat v5 style)", {
  markers <- suppressWarnings(cosg(
    object = pbmc_small,
    groups = 'all',
    assay = 'RNA',
    layer = 'data',
    mu = 1,
    n_genes_user = 100
  ))

  # Check output structure

expect_true(is.list(markers))
  expect_true("names" %in% names(markers))
  expect_true("scores" %in% names(markers))

  # Check dimensions
  expect_equal(nrow(markers$names), 100)
  expect_equal(nrow(markers$scores), 100)
  expect_equal(ncol(markers$names), 3)  # 3 clusters in pbmc_small
  expect_equal(ncol(markers$scores), 3)

  # Check that we get gene names and numeric scores
  expect_true(is.character(markers$names[1, 1]))
  expect_true(is.numeric(markers$scores[1, 1]))
})

test_that("COSG works with slot parameter (Seurat v3/v4 backward compatibility)", {
  markers <- suppressWarnings(cosg(
    object = pbmc_small,
    groups = 'all',
    assay = 'RNA',
    slot = 'data',
    mu = 1,
    n_genes_user = 100
  ))

  # Check output structure
  expect_true(is.list(markers))
  expect_equal(nrow(markers$names), 100)
  expect_equal(ncol(markers$names), 3)
})

test_that("layer parameter takes precedence over slot parameter", {
  # When both are provided, layer should be used
  markers_layer <- suppressWarnings(cosg(
    object = pbmc_small,
    groups = 'all',
    assay = 'RNA',
    layer = 'data',
    mu = 1,
    n_genes_user = 50
  ))

  markers_both <- suppressWarnings(cosg(
    object = pbmc_small,
    groups = 'all',
    assay = 'RNA',
    slot = 'counts',  # Different slot
    layer = 'data',   # But layer takes precedence
    mu = 1,
    n_genes_user = 50
  ))

  # Results should be identical when layer is the same
  expect_equal(markers_layer$names, markers_both$names)
})

test_that("COSG works with selected groups", {
  markers <- suppressWarnings(cosg(
    object = pbmc_small,
    groups = c('0', '2'),
    assay = 'RNA',
    layer = 'data',
    mu = 1,
    n_genes_user = 50
  ))

  # Should only have 2 columns for the 2 selected groups
  expect_equal(ncol(markers$names), 2)
  expect_equal(ncol(markers$scores), 2)
  expect_true(all(colnames(markers$names) %in% c('0', '2')))
})

test_that("COSG works with different mu values", {
  markers_mu1 <- suppressWarnings(cosg(
    object = pbmc_small,
    layer = 'data',
    mu = 1,
    n_genes_user = 50
  ))

  markers_mu10 <- suppressWarnings(cosg(
    object = pbmc_small,
    layer = 'data',
    mu = 10,
    n_genes_user = 50
  ))

  # Different mu values should produce different results
  expect_false(identical(markers_mu1$scores, markers_mu10$scores))
})

test_that("COSG works with remove_lowly_expressed parameter", {
  markers <- suppressWarnings(cosg(
    object = pbmc_small,
    layer = 'data',
    mu = 1,
    n_genes_user = 50,
    remove_lowly_expressed = TRUE,
    expressed_pct = 0.1
  ))

  expect_equal(nrow(markers$names), 50)
})

test_that("COSG errors on single cluster", {
  expect_error(
    cosg(object = pbmc_small, groups = c('0')),
    "Cannot perform marker gene identification on a single cluster"
  )
})
