context("Identification of Marker Genes by COSG")

library('COSG')


markers <- suppressWarnings(cosg(object = pbmc_small, alpha = 0.1,n_genes_user=100))


test_that("COSG works as expected", {
  expect_error(cosg(object = pbmc_small), NA)
  expect_error(cosg(object = pbmc_small, alpha = 0.1), NA)
  expect_error(cosg(object = pbmc_small, alpha = 0.1,n_genes_user=100), NA)
    
  expect_equal(dim(markers$names), c(100,3))
  expect_equal(dim(markers$scores), c(100,3))

  expect_equal(markers$names[1,1], 'CCL5')
  expect_equal(markers$scores[1,1], 0.687876833896255,tolerance = 1e-6)
})
#> Test passed 🎉