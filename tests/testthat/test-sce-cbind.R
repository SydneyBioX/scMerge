context("Testing sce_cbind function")

## Loading required data
library(SingleCellExperiment)
data("example_sce", package = "scMerge")

## A simple case: The sce_mESC example data contains 5 different batches.  We will split this data by batches first and combine it
batch_names = unique(example_sce$batch)
sce_split = list(example_sce[, example_sce$batch == "batch2"], example_sce[, example_sce$batch == "batch3"])

## Testing if the combination yields the same outputs
sce_combine = sce_cbind(sce_list = sce_split, 
                        batch_names = batch_names, 
                        cut_off_batch = 0, 
                        cut_off_overall = 0, 
                        method = "union", 
                        exprs = c("counts", "logcounts"), 
                        colData_names = c("cellTypes"))

expect_identical(dim(example_sce), dim(sce_combine))
expect_identical(colData(example_sce), colData(sce_combine))
expect_identical(colData(example_sce)$batch, colData(sce_combine)$batch)


## Testing if the combination yields the same outputs
sce_combine = sce_cbind(sce_list = sce_split, 
                        batch_names = batch_names, 
                        cut_off_batch = 0, 
                        cut_off_overall = 0, 
                        method = "intersect", 
                        exprs = c("counts", "logcounts"), 
                        colData_names = c("cellTypes"))

## Mis-specified method argument
expect_error(sce_cbind(sce_list = sce_split, method = "abcdefg"))


################
library(SummarizedExperiment)

da_sce_split = sce_split
assay(da_sce_split[[1]], "counts") = DelayedArray::DelayedArray(assay(sce_split[[1]], "counts"))
assay(da_sce_split[[1]], "logcounts") = DelayedArray::DelayedArray(assay(sce_split[[1]], "logcounts"))

assay(da_sce_split[[2]], "counts") = DelayedArray::DelayedArray(assay(sce_split[[2]], "counts"))
assay(da_sce_split[[2]], "logcounts") = DelayedArray::DelayedArray(assay(sce_split[[2]], "logcounts"))

sce_cbind(sce_list = da_sce_split, 
          batch_names = batch_names, 
          cut_off_batch = 0, 
          cut_off_overall = 0, 
          method = "intersect", 
          exprs = c("counts", "logcounts"), 
          colData_names = c("cellTypes"))

expect_error(
  sce_cbind(sce_list = da_sce_split, 
            batch_names = batch_names, 
            cut_off_batch = 0, 
            cut_off_overall = 0, 
            method = "union", 
            exprs = c("counts", "logcounts"), 
            colData_names = c("cellTypes")))
