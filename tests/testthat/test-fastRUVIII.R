context("fastRUVIII")

## Generate a simulated data to test on the speed of our improved algorithm
L = scMerge::ruvSimulate(m = 800, n = 1000, nc = 50, nRep = 10)
Y = L$Y; M = L$M; ctl = L$ctl

t1 = Sys.time()
improved1 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, fast_svd = FALSE)
t2 = Sys.time()
improved2 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, fast_svd = TRUE, rsvd_prop = 0.1)
t3 = Sys.time()
old = ruv::RUVIII(Y = Y, M = M, ctl = ctl, k = 20)
t4 = Sys.time()

## improved1 used C++ matrix operations, hence it should be faster than native R matrix operations
expect_equal(improved1, old)
expect_lt(as.numeric(t2 - t1, units = "secs"),
          as.numeric(t4 - t3, units = "secs"))

## improved2 used RSVD, hence it should be faster than old
expect_lt(as.numeric(t3 - t2, units = "secs"),
          as.numeric(t4 - t3, units = "secs"))

## improved2 used RSVD, hence it should be faster than improved1
expect_lt(as.numeric(t3 - t2, units = "secs"),
          as.numeric(t2 - t1, units = "secs"))
