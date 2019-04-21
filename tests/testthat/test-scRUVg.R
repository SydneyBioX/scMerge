context("Test scRUVg")
set.seed(1234)
L = scMerge::ruvSimulate(m = 80, n = 1000, nc = 50, nCelltypes = 10)
Y = L$Y; ctl = L$ctl; X = L$X
ruvgRes = scMerge::scRUVg(Y = Y, ctl = ctl, k = 20)
old = ruv::RUV2(Y = Y, X = X, ctl = ctl, k = 20)
expect_equal(abs(ruvgRes$W), abs(old$W))
