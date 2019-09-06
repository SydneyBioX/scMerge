context("Handling sparse matrix in sce")
library(HDF5Array)

set.seed(12345)
L = ruvSimulate(m = 1000, n = 200, nc = 50, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

# pryr::object_size(L)

system.time({
  non_delayed = scMerge(
    sce_combine = L,
    ctl = paste0("gene",1:100),
    exprs = "logcounts",
    hvg_exprs = "counts",
    kmeansK = c(3, 3),
    cell_type = L$cellTypes,
    replicate_prop = 1,
    assay_name = 'scMerge')
})

# assay(L, "counts") = DelayedArray::DelayedArray(seed = assay(L, "counts"))
# assay(L, "logcounts") = DelayedArray::DelayedArray(seed = assay(L, "logcounts"))

assay(L, "counts") = as(assay(L, "counts"), "HDF5Matrix")
assay(L, "logcounts") = as(assay(L, "logcounts"), "HDF5Matrix")
# pryr::object_size(L)

system.time({
  delayed = scMerge(
    sce_combine = L,
    ctl = paste0("gene",1:100),
    exprs = "logcounts", 
    hvg_exprs = "counts",
    cell_type = L$cellTypes,
    replicate_prop = 1,
    assay_name = 'scMerge')
})

# pryr::object_size(delayed)
