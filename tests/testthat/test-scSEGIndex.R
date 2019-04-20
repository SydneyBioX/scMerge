context("Test scSEGIndex")
set.seed(1)
L = ruvSimulate(m = 50, n = 150, nc = 10, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = TRUE)
scSEGIndex(exprsMat = SummarizedExperiment::assay(L, 'counts'))
