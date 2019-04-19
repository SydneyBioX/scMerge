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

## Simulated data
set.seed(12345)
L = ruvSimulate(m = 100, n = 1000, nc = 400, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

sce_mESC <- scMerge(
  sce_combine = L,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = L$cellTypes,
  replicate_prop = 1,
  assay_name = 'scMerge')