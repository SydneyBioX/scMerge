context("Testing sce_cbind function")

## Loading required data
require(SingleCellExperiment)
data("sce_mESC", package = "scMerge.data")

## A simple case: The sce_mESC example data contains 5 different batches. We will split this data by batches first and combine it
batch_names = unique(sce_mESC$batch)
sce_split = list(
  sce_mESC[,sce_mESC$batch == batch_names[1]],
  sce_mESC[,sce_mESC$batch == batch_names[2]],
  sce_mESC[,sce_mESC$batch == batch_names[3]],
  sce_mESC[,sce_mESC$batch == batch_names[4]],
  sce_mESC[,sce_mESC$batch == batch_names[5]]
)

## Testing if the combination yields the same outputs
sce_combine = sce_cbind(sce_list = sce_split,
                        batch_names = batch_names,
                        cut_off_batch = 0,
                        cut_off_overall = 0,
                        method = NULL,
                        exprs = c("counts", "logcounts"),
                        colData_names = c("cellTypes"))
# expect_identical(sce_mESC, sce_combine)
# expect_equal(sce_mESC@assays, sce_combine@assays)
# expect_identical(sce_mESC@int_metadata, sce_combine@int_metadata)
# expect_identical(sce_mESC@reducedDims, sce_combine@reducedDims)
expect_identical(dim(sce_mESC), dim(sce_combine))
expect_identical(colData(sce_mESC), colData(sce_combine))
expect_identical(colData(sce_mESC)$batch, colData(sce_combine)$batch)

expect_error(
  sce_cbind(sce_list = sce_split,
            method = "all")
)
