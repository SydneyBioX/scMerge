# library(microbenchmark)
# library(DelayedArray)
# library(Matrix)
# library(pryr)
# library(HDF5Array)
# setAutoBPPARAM(SerialParam())
# DelayedArray:::set_verbose_block_processing(TRUE)
# 
# set.seed(1)
# L = ruvSimulate(m = 1000, n = 20000, nCelltypes = 10, lambda = 0.1)
# Y = L$Y
# M = L$M
# sY = as(Y, "dgeMatrix")
# sM = as(M, "dgeMatrix")
# dY = DelayedArray(Y)
# dM = DelayedArray(M)
# fY = as(Y, "HDF5Array")
# fM = as(M, "HDF5Array")
# 
# pryr::object_size(Y)
# pryr::object_size(sY)
# pryr::object_size(dY)
# pryr::object_size(fY)
# 
# 
# microbenchmark::microbenchmark(
#   scMerge::eigenResidop(Y, M), 
#   my_residop(Y, M),
#   ruv::residop(Y, M),
#   my_residop(sY, sM),
#   my_residop(dY, dM),
#   my_residop(fY, fM),
#   times = 10)
