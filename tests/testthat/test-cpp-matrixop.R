context("Testing speed of cpp matrix operations")

## Generate a simulated data to test on the speed of our improved matrix
## operations
L = ruvSimulate(m = 100, n = 2000, nRep = 10, lambda = 0.1)
Y = L$Y
M = L$M


t1 = Sys.time()
improved1 = unname(scMerge::eigenResidop(Y, M))
t2 = Sys.time()
old = unname(ruv::residop(Y, M))
t3 = Sys.time()

## improved1 used C++ matrix operations, hence it should be faster than
## native R matrix operations
expect_equal(improved1, old, tolerance = 1e-07)
expect_lt(as.numeric(t2 - t1, units = "secs"), as.numeric(t3 - t2, units = "secs"))
