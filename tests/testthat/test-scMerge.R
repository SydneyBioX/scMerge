context("Test scMerge")

data('example_sce', package = 'scMerge')
## Previously computed stably expressed genes
data('segList_ensemblGeneID', package = 'scMerge')
## Running an example data with minimal inputs
sce_mESC <- scMerge(
  sce_combine = example_sce,
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = 'scMerge')

## Simulated data, testing on sce = TRUE option and supervised scMerge
set.seed(12345)
L = ruvSimulate(m = 100, n = 1000, nc = 400, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

L_result1 <- scMerge(
  sce_combine = L,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = L$cellTypes,
  replicate_prop = 1,
  assay_name = 'scMerge')


ruvK = c(10, 20)
L_result2 <- scMerge(
  sce_combine = L,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = L$cellTypes,
  replicate_prop = 1,
  # return_all_RUV = TRUE, 
  ruvK = ruvK,
  assay_name = paste0("ruv", ruvK))


########## Checking error handling

## assay_name NULL (by default it is NULL)

expect_error(
  sce_mESC <- scMerge(
    sce_combine = example_sce,
    ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
    kmeansK = c(3, 3))
)


## duplicated cell names
L2 = L
colnames(L2)[1:2] = "cellDups"

expect_error(
  scMerge(
    sce_combine = L2,
    ctl = paste0("gene",1:100),
    kmeansK = c(3, 3),
    cell_type = L$cellTypes,
    replicate_prop = 1,
    assay_name = 'scMerge')
)

## If return_all_RUV = TRUE, then length(ruvK) must match length(assay_name)

expect_error(
  scMerge(
    sce_combine = example_sce,
    ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
    kmeansK = c(3, 3),
    assay_name = 'scMerge', 
    return_all_RUV = TRUE, 
    ruvK = c(10, 20))
)


## Both exprs matrix and hvg_exprs matrix must be base::matrix objects 
## (This will be updated in later versions of scMerge)
L3 = L
SingleCellExperiment::counts(L3) = 
  Matrix::Matrix(SingleCellExperiment::counts(L3), sparse = TRUE)

expect_error(
  scMerge(
    sce_combine = L3,
    ctl = paste0("gene",1:100),
    assay_name = 'scMerge')
)

L4 = L
SingleCellExperiment::logcounts(L4) = 
  Matrix::Matrix(SingleCellExperiment::logcounts(L4), sparse = TRUE)

expect_error(
  scMerge(
    sce_combine = L4,
    ctl = paste0("gene",1:100),
    assay_name = 'scMerge')
)

## Checking batch
## Making batch as characters
example_sce2 = example_sce
example_sce2$batch = as.character(example_sce2$batch)
scMerge(
  sce_combine = example_sce2,
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = 'scMerge')

## Inputting NULL batch
example_sce2$batch = NULL

expect_error(
  scMerge(
    sce_combine = example_sce2,
    ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
    kmeansK = c(3, 3),
    assay_name = 'scMerge')
)
