context("Testing speed of cpp matrix operations")

## Generate a simulated data to test on the speed of our improved matrix operations
set.seed(1)
L = ruvSimulate(m = 100, n = 2000, nCelltypes = 10, lambda = 0.1)
Y = L$Y
M = L$M


improved1 = unname(scMerge::eigenResidop(Y, M))
improved2 = unname(scMerge::eigenMatMult(t(M), Y))
old1 = unname(ruv::residop(Y, M))
old2 = unname(t(M) %*% Y)


## improved1 and 2 used C++ matrix operations, hence it should be faster than native R matrix operations
expect_equal(improved1, old1, tolerance = 1e-06)
expect_equal(improved2, old2, tolerance = 1e-06)
