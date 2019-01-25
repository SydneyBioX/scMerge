context("Testing if parallelisation produce identical results")

## Loading required data
require(SingleCellExperiment)
data("sce_mESC", package = "scMerge.data")
data("segList_ensemblGeneID")

t1 = Sys.time()
pF <- scMerge(sce_mESC,
              ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
              kmeansK = c(1,3,3,1,1),
              assay_name = "scMerge_classical",
              parallel = FALSE)
t2 = Sys.time()
t2 - t1

t3 = Sys.time()
pT <- scMerge(sce_mESC,
              ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
              kmeansK = c(1,3,3,1,1),
              assay_name = "scMerge_classical",
              parallel = TRUE)
t4 = Sys.time()
t4 - t3

expect_identical(
  assay(pF, "scMerge_classical"),
  assay(pT, "scMerge_classical")
)

all.equal(pT, pF)
