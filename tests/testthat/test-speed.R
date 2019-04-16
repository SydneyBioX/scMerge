context("Test speed")
set.seed(12345)
L = ruvSimulate(m = 1000, n = 10000, nc = 400, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
# Y = L$Y; M = L$M; ctl = L$ctl;
# dim(Y)
# max(dim(Y))/min(dim(Y))
# system.time(svd(Y))
# system.time(svd(eigenMatMult(Y, t(Y))))
# 
# system.time(scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, fast_svd = FALSE))
# system.time(scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, fast_svd = TRUE))
data("example_sce", package = "scMerge")
data("segList_ensemblGeneID", package = "scMerge")
dim(example_sce)
system.time(scMerge(
  sce_combine = example_sce,
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = 'scMerge'))