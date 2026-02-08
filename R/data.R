#' A small example Seurat object for testing
#'
#' A subsetted version of a 10X Genomics PBMC dataset, containing 80 cells
#' and 230 genes. This is the same dataset as SeuratObject::pbmc_small.
#'
#' @format A Seurat object with:
#' \describe{
#'   \item{assays}{RNA assay with counts and data layers}
#'   \item{meta.data}{Cell metadata including cell type annotations}
#' }
#' @source SeuratObject package
#' @examples
#' data(pbmc_small)
#' pbmc_small
"pbmc_small"
