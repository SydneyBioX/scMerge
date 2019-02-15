context("Testing identical results when permuting columns")

library(SingleCellExperiment)
library(scater)
library(scMerge)

## Loading example data
data("example_sce", package = "scMerge")
## Previously computed stably expressed genes
data("segList_ensemblGeneID", package = "scMerge")
## Running an example data with minimal inputs
example_sce1 = example_sce
res1 = scMerge(sce_combine = example_sce1, ctl = segList_ensemblGeneID$mouse$mouse_scSEG, 
    kmeansK = c(3, 3), assay_name = "scMerge1", replicate_prop = 1)
###################################### 
set.seed(2)
index = sample(seq_len(ncol(example_sce)))
example_sce2 = example_sce[, index]
res2 = scMerge(sce_combine = example_sce2, ctl = segList_ensemblGeneID$mouse$mouse_scSEG, 
    kmeansK = c(3, 3), assay_name = "scMerge2", replicate_prop = 1)

# scater::plotPCA(res1, colour_by = 'cellTypes', shape = 'batch', run_args =
# list(exprs_values = 'scMerge1'), add_ticks = FALSE)

# scater::plotPCA(res2, colour_by = 'cellTypes', shape = 'batch', run_args =
# list(exprs_values = 'scMerge2'), add_ticks = FALSE)


## Checking if the ourput matrices are numerically equal
expect_equal(assay(res1, "scMerge1"), assay(res2, "scMerge2")[, order(index)])


## Checking if the replication matrices are numerically equal up to a
## permutation
repMat1 = unname(metadata(res1)$scRep_res)
repMat2 = unname(metadata(res2)$scRep_res[order(index), ])

expect_true(identical(repMat1[, 1], repMat2[, 1]) | identical(repMat1[, 1], 
    repMat2[, 2]) | identical(repMat1[, 1], repMat2[, 3]))

expect_true(identical(repMat1[, 2], repMat2[, 1]) | identical(repMat1[, 2], 
    repMat2[, 2]) | identical(repMat1[, 2], repMat2[, 3]))

expect_true(identical(repMat1[, 3], repMat2[, 1]) | identical(repMat1[, 3], 
    repMat2[, 2]) | identical(repMat1[, 3], repMat2[, 3]))
