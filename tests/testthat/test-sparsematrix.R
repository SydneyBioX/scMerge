context("Test sparseArray")

library(Matrix)

set.seed(12345)
L = ruvSimulate(m = 500, n = 20000, nc = 400, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = TRUE)

sce_matrix <- scMerge(
  sce_combine = L,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = L$cellTypes,
  replicate_prop = 1,
  assay_name = 'matrix_output')

counts = assay(sce_matrix, "counts")
logcounts = assay(sce_matrix, "logcounts")
################################################

sce_sp = sce_matrix
assay(sce_sp, "counts") = as(counts, "dgCMatrix")
assay(sce_sp, "logcounts") = as(logcounts, "dgCMatrix")

sce_sp <- scMerge(
  sce_combine = sce_sp,
  ctl = paste0("gene",1:100),
  kmeansK = c(3, 3),
  cell_type = sce_sp$cellTypes,
  replicate_prop = 1,
  assay_name = 'sp_output')

expect_identical(as.matrix(assay(sce_sp, "sp_output")), 
                 assay(sce_matrix, "matrix_output"))

#################################################################

L = ruvSimulate(m = 1000, n = 20000, nc = 400, 
                nCelltypes = 3, nBatch = 2, 
                lambda = 0.1, sce = FALSE)
Y = L$Y
M = L$M
Y_sp = as(Y, "dgCMatrix")
M_sp = as(M, "dgCMatrix")
mean(Y == 0)
mean(M == 0)

sp_residop = function (A, B)
{
  return(A - B %*% base::solve(t(B) %*% B) %*% t(B) %*% A)
}

# microbenchmark::microbenchmark(
#   t(Y) %*% M,
#   eigenMatMult(t(Y), M),
#   t(Y_sp) %*% M_sp,
#   eigenSpMatMult(t(Y_sp), M_sp),
#   times = 20)

# microbenchmark::microbenchmark(
#   ruv::residop(Y, M),
#   eigenResidop(Y, M),
#   eigenSpResidop(Y_sp, M_sp),
#   sp_residop(Y_sp, M_sp),
#   times = 5)

profvis::profvis({
  ruv::residop(Y, M)
  eigenResidop(Y, M)
  sp_residop(Y_sp, M_sp)
  eigenSpResidop(Y_sp, M_sp)
})


profvis::profvis({
  t(Y) %*% M
  eigenMatMult(t(Y), M)
  t(Y_sp) %*% M_sp
  eigenSpMatMult(t(Y_sp), M_sp)
})