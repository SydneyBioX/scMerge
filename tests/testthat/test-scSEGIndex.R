## Loading example data
data('example_sce', package = 'scMerge')
## subsetting genes to illustrate usage.
exprs_mat = SummarizedExperiment::assay(example_sce, 'logcounts')[1:110, 1:20]
set.seed(1)
scSEGIndex(exprs_mat = exprs_mat)
