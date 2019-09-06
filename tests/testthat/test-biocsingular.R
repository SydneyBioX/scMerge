context("Test biocsingular")

library(HDF5Array)
data('example_sce', package = 'scMerge')
## Previously computed stably expressed genes
data('segList_ensemblGeneID', package = 'scMerge')
## Running an example data with minimal inputs
sce_mESC <- scMerge(
  sce_combine = example_sce,
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = 'scMerge')

example_sce_da = example_sce

# stats::cor(matrix(rnorm(100), 10, 10) %>% as(., "HDF5Array") %>% as.matrix, 
#            matrix(rnorm(100), 10, 10) %>% as(., "HDF5Array") %>% as.matrix)

assay(example_sce_da, "counts") = as(assay(example_sce_da, "counts"), "HDF5Array")
assay(example_sce_da, "logcounts") = as(assay(example_sce_da, "logcounts"), "HDF5Array")

example_sce_da <- scMerge(
  sce_combine = example_sce_da,
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = 'scMerge')
