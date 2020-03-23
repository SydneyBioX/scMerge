context("Test DelayedArray")

library(DelayedArray)
library(BiocParallel)
BiocParallel::register(BPPARAM = BiocParallel::SerialParam())

set.seed(12345)
L = ruvSimulate(m = 100, n = 1000, nc = 100, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

sce_matrix <- scMerge(
  sce_combine = L,
  ctl = paste0("gene",1:10),
  kmeansK = c(3, 3),
  cell_type = L$cellTypes,
  assay_name = 'matrix_output')

counts = assay(L, "counts")
logcounts = assay(L, "logcounts")
################################################

sce_hdf = L
assay(sce_hdf, "counts") = DelayedArray(counts)
assay(sce_hdf, "logcounts") = DelayedArray(logcounts)
# setRealizationBackend("HDF5Array")
# getRealizationBackend()
# library(BiocParallel)
# DelayedArray::setAutoBPPARAM(BPPARAM = BiocParallel::MulticoreParam(workers = 5))
# getAutoBPPARAM()

sce_hdf <- scMerge(
    sce_combine = sce_hdf,
    ctl = paste0("gene",1:10),
    kmeansK = c(3, 3),
    cell_type = sce_hdf$cellTypes,
    assay_name = 'hdf_output')

sce_hdf

expect_equal(as.matrix(assay(sce_hdf, "hdf_output")), 
                 assay(sce_matrix, "matrix_output"))