context("Test parallelisation")

# library(doParallel)
# library(BiocParallel)
# register(bpstart(MulticoreParam(workers=5)))
# profvis::profvis({
#   BiocSingular::runRandomSVD(x = Y0, k = 10, BPPARAM = SerialParam())
#   BiocSingular::runRandomSVD(x = Y0, k = 10, BPPARAM = BiocParallel::MulticoreParam(workers = 5))
# })
