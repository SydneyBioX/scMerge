context("Testing on scReplicate")

library(SingleCellExperiment)
library(scMerge)
library(scMerge.data)

## Loading example data
data("sce_mESC", package = "scMerge.data")

#################### Testing if the fast_svd option yields identical result ####################
set.seed(1)
t1 <- Sys.time()
scRep_fastF <- scReplicate(
  sce = sce_mESC,
  batch = sce_mESC$batch,
  kmeansK = c(1,3,3,1,1),
  fast_svd = FALSE)
t2 <- Sys.time()


set.seed(1)
scRep_fastT <- scReplicate(
  sce = sce_mESC,
  batch = sce_mESC$batch,
  kmeansK = c(1,3,3,1,1),
  fast_svd = TRUE)
t3 <- Sys.time()

expect_identical(scRep_fastF, scRep_fastT)
## The time for fast_svd = FALSE should be greater than when fast_svd = TRUE
expect_gt(as.numeric(t2 - t1, units = "secs"),
          as.numeric(t3 - t2, units = "secs"))



#################### Testing if the parallel option yields identical result ####################
## mESC data is not large enough such that the parallelisation provides an improvement in time. 
## Ask YXL to test this code on an alterntive data.
# t1 <- Sys.time()
# scRep_pT <- scReplicate(
#   sce = sce_mESC,
#   batch = sce_mESC$batch,
#   kmeansK = c(1, 3, 3, 1, 1),
#   fast_svd = TRUE,
#   parallelParam = SnowParam(workers = 5, type = "FORK"))
# t2 <- Sys.time()
# 
# 
# t3 <- Sys.time()
# scRep_pF <- scReplicate(
#   sce = sce_mESC,
#   batch = sce_mESC$batch,
#   kmeansK = c(1, 3, 3, 1, 1),
#   fast_svd = TRUE,
#   parallelParam = BiocParallel::SerialParam())
# t4 <- Sys.time()
# 
# expect_identical(scRep_pF, scRep_pT)
# 
# as.numeric(t2 - t1, units = "secs")
# as.numeric(t4 - t3, units = "secs")
# # expect_gt(as.numeric(t2 - t1, units = "secs"),
# #           as.numeric(t3 - t2, units = "secs"))
