context("Test on scRUVIII")

set.seed(1)
L = scMerge::ruvSimulate(m = 200, n = 1000, nc = 50, nRep = 10)
Y = L$Y; M = L$M; ctl = L$ctl; batch = L$dataSource;
res = scRUVIII(Y = Y, M = M, ctl = ctl, k = 20, batch = batch)

## When multiple k's are supplied, we have a message indicating a selection of k using silhoutte coefficients. 
expect_message(
  scRUVIII(Y = Y, M = M, ctl = ctl, k = c(10, 20), batch = batch)
)
