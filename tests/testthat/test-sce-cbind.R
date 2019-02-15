context("Testing sce_cbind function")

## Loading required data
library(SingleCellExperiment)
data("example_sce", package = "scMerge")

## A simple case: The sce_mESC example data contains 5 different batches.  We
## will split this data by batches first and combine it
batch_names = unique(example_sce$batch)
sce_split = list(example_sce[, example_sce$batch == "batch2"], example_sce[, 
    example_sce$batch == "batch3"])

## Testing if the combination yields the same outputs
sce_combine = sce_cbind(sce_list = sce_split, batch_names = batch_names, cut_off_batch = 0, 
    cut_off_overall = 0, method = NULL, exprs = c("counts", "logcounts"), colData_names = c("cellTypes"))
# expect_identical(example_sce, sce_combine)
# expect_equal(example_sce@assays, sce_combine@assays)
# expect_identical(example_sce@int_metadata, sce_combine@int_metadata)
# expect_identical(example_sce@reducedDims, sce_combine@reducedDims)
expect_identical(dim(example_sce), dim(sce_combine))
expect_identical(colData(example_sce), colData(sce_combine))
expect_identical(colData(example_sce)$batch, colData(sce_combine)$batch)

## Mis-specified method argument
expect_error(sce_cbind(sce_list = sce_split, method = "all"))
