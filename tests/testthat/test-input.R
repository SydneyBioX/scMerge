context("Inputs of scMerge")

## Input must be an SingleCellExperiment object
expect_error(
  scMerge(
    sce_combine = matrix(c(0,1,1,1))
    )
  )

## Input sce must have unique cell names
data("example_sce", package = "scMerge")
example_sce_colnames = example_sce
colnames(example_sce_colnames) = rep("Cell", ncol(example_sce_colnames))
expect_error(
  scMerge(
    sce_combine = example_sce_colnames
  )
)

## Must have a value for "exprs" argument, the default is "logcounts"
expect_error(
  scMerge(
    sce_combine = example_sce,
    exprs = NULL
  )
)


## Must have a value for "assay_name" argument
expect_error(
  scMerge(
    sce_combine = example_sce
  )
)


## Must have valid ctl values. 
expect_error(
  scMerge(
    sce_combine = example_sce,
    assay_name = "scmerge"
  )
)

expect_error(
  scMerge(
    sce_combine = example_sce,
    assay_name = "scmerge",
    ctl = c("apple", "banana")
  )
)
