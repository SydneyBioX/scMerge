context("Test biocsingular")

library(HDF5Array)

set.seed(12345)
L = ruvSimulate(m = 100, n = 1000, nc = 400, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

sce_matrix <- scMerge(
  sce_combine = L,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = L$cellTypes,
  replicate_prop = 1,
  assay_name = 'matrix_output')

counts = assay(sce_matrix, "counts")
logcounts = assay(sce_matrix, "logcounts")
sce_da = sce_matrix
assay(sce_da, "counts") = as(counts, "HDF5Array")
assay(sce_da, "logcounts") = as(logcounts, "HDF5Array")

sce_da <- scMerge(
  sce_combine = sce_da,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = sce_da$cellTypes,
  replicate_prop = 1,
  assay_name = 'da_output')

expect_identical(as.matrix(assay(sce_da, "da_output")), 
                 assay(sce_matrix, "matrix_output"))



sce_Matrix = sce_matrix
assay(sce_Matrix, "counts") = DelayedArray(Matrix::Matrix(counts))
assay(sce_Matrix, "logcounts") = DelayedArray(Matrix::Matrix(logcounts))

sce_Matrix <- scMerge(
  sce_combine = sce_Matrix,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = sce_Matrix$cellTypes,
  replicate_prop = 1,
  assay_name = 'dgC_output')

expect_identical(as.matrix(assay(sce_Matrix, "dgC_output")), 
          assay(sce_matrix, "matrix_output"))
