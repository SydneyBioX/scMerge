context("Test scRUVIII")
set.seed(1234)
L = ruvSimulate(m = 200, n = 1000, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
Y = t(log2(L$Y + 1L)); M = L$M; ctl = L$ctl; batch = L$batch;
res = scRUVIII(Y = Y, M = M, ctl = ctl, k = c(5, 10, 15, 20), batch = batch)
