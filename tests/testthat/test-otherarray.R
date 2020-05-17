context("Test DelayedArray and SparseArray")

library(Matrix)
library(DelayedArray)
library(BiocParallel)
library(HDF5Array)
# library(BiocParallel)
# DelayedArray::setAutoBPPARAM(BPPARAM = BiocParallel::MulticoreParam(workers = 5))
# DelayedArray:::set_verbose_block_processing(TRUE)
# BiocParallel::register(BPPARAM = BiocParallel::SerialParam())
# setAutoBPPARAM(SerialParam())

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
################ Sparse ################
sce_sp = L
assay(sce_sp, "counts") = as(counts, "dgeMatrix")
assay(sce_sp, "logcounts") = as(logcounts, "dgeMatrix")
expect_equal(counts, as.matrix(assay(sce_sp, "counts")))
expect_equal(logcounts, as.matrix(assay(sce_sp, "logcounts")))

sce_sp <- scMerge(
  sce_combine = sce_sp,
  ctl = paste0("gene",1:10),
  kmeansK = c(3, 3),
  cell_type = sce_sp$cellTypes,
  assay_name = 'sp_output')

expect_equal(as.matrix(assay(sce_sp, "sp_output")), 
             assay(sce_matrix, "matrix_output"))
################ DelayedArray ################
sce_da = L
assay(sce_da, "counts") = DelayedArray(counts)
assay(sce_da, "logcounts") = DelayedArray(logcounts)

sce_da <- scMerge(
  sce_combine = sce_da,
  ctl = paste0("gene",1:10),
  kmeansK = c(3, 3),
  cell_type = sce_da$cellTypes,
  assay_name = 'da_output')

sce_da

expect_equal(as.matrix(assay(sce_da, "da_output")), 
             assay(sce_matrix, "matrix_output"))
################ HDF5Array ################
sce_hdf = L
assay(sce_hdf, "counts") = as(counts, "HDF5Array")
assay(sce_hdf, "logcounts") = as(logcounts, "HDF5Array")

sce_hdf <- scMerge(
  sce_combine = sce_hdf,
  ctl = paste0("gene",1:10),
  kmeansK = c(3, 3),
  cell_type = sce_hdf$cellTypes,
  assay_name = 'hdf_output')

sce_hdf

expect_equal(as.matrix(assay(sce_hdf, "hdf_output")), 
             assay(sce_matrix, "matrix_output"))