context("Test fastRUVIII")
library(BiocSingular)
set.seed(12345)
L = ruvSimulate(m = 400, n = 500, nc = 400, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
Y = L$Y
M = L$M
ctl = L$ctl

old = ruv::RUVIII(Y = Y, M = M, ctl = ctl, k = 20)

improved1 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, BSPARAM = ExactParam())

improved2 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, BSPARAM = RandomParam(), svd_prop = 0.3)

expect_equal(old, improved1)
expect_equal(improved1, improved2, tol = 0.01)
# expect_lt(as.numeric(t3 - t2, units = 'secs'), as.numeric(t2 - t1, units = 'secs')) expect_lt(as.numeric(t4 - t3, units = 'secs'),
# as.numeric(t3 - t2, units = 'secs'))
