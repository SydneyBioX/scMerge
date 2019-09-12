context("Test sparseArray")

library(Matrix)

set.seed(12345)
L = ruvSimulate(m = 100, n = 2000, nc = 400, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

# profvis::profvis(
sce_matrix <- scMerge(
  sce_combine = L,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = L$cellTypes,
  replicate_prop = 1,
  assay_name = 'matrix_output')
# )
counts = assay(L, "counts")
logcounts = assay(L, "logcounts")
################################################

sce_sp = L
assay(sce_sp, "counts") = as(counts, "dgCMatrix")
assay(sce_sp, "logcounts") = as(logcounts, "dgCMatrix")

# profvis::profvis(
sce_sp <- scMerge(
  sce_combine = sce_sp,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = sce_sp$cellTypes,
  replicate_prop = 1,
  assay_name = 'sp_output')
# )

expect_equal(as.matrix(assay(sce_sp, "sp_output")), 
                 assay(sce_matrix, "matrix_output"))

####################################################################
sce_sp_da = L
assay(sce_sp_da, "counts") = DelayedArray::DelayedArray(as(counts, "dgCMatrix"))
assay(sce_sp_da, "logcounts") = DelayedArray::DelayedArray(as(logcounts, "dgCMatrix"))

# profvis::profvis(
sce_sp_da <- scMerge(
  sce_combine = sce_sp_da,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = sce_sp_da$cellTypes,
  replicate_prop = 1,
  assay_name = 'sp_da_output')
# )

expect_equal(as.matrix(assay(sce_sp_da, "sp_da_output")),
             assay(sce_matrix, "matrix_output"))