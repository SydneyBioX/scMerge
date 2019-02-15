context("Testing if parallelisation produce identical results")

## Loading required data
require(SingleCellExperiment)
data("example_sce", package = "scMerge")
data("segList_ensemblGeneID", package = "scMerge")

t1 = Sys.time()
pF <- scMerge(example_sce, ctl = segList_ensemblGeneID$mouse$mouse_scSEG, kmeansK = c(3, 
    3), assay_name = "scMerge_classical", parallel = FALSE)
t2 = Sys.time()
t2 - t1


t3 = Sys.time()
pT <- scMerge(example_sce, ctl = segList_ensemblGeneID$mouse$mouse_scSEG, kmeansK = c(3, 
    3), assay_name = "scMerge_classical", parallel = TRUE, parallelParam = SnowParam(workers = 8, 
    type = "FORK"))
t4 = Sys.time()
t4 - t3

expect_identical(assay(pF, "scMerge_classical"), assay(pT, "scMerge_classical"))

all.equal(pT, pF)
