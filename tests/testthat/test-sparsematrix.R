context("Test sparseArray")

library(Matrix)

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
################################################

sce_sp = sce_matrix
assay(sce_sp, "counts") = as(counts, "dgCMatrix")
assay(sce_sp, "logcounts") = as(logcounts, "dgCMatrix")

sce_sp <- scMerge(
  sce_combine = sce_sp,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = sce_sp$cellTypes,
  replicate_prop = 1,
  assay_name = 'sp_output')

expect_identical(as.matrix(assay(sce_sp, "sp_output")), 
                 assay(sce_matrix, "matrix_output"))
