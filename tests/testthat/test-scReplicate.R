context("Testing on scReplicate")

library(SingleCellExperiment)
library(scMerge)

## Loading example data
data("example_sce", package = "scMerge")

#################### Testing if the fast_svd option yields identical result
#################### ####################
set.seed(1)
t1 <- Sys.time()
scRep_fastF <- scReplicate(sce = example_sce, batch = example_sce$batch, kmeansK = c(3, 
    3), fast_svd = FALSE)
t2 <- Sys.time()


set.seed(1)
scRep_fastT <- scReplicate(sce = example_sce, batch = example_sce$batch, kmeansK = c(3, 
    3), fast_svd = TRUE)
t3 <- Sys.time()

## Since the clustering could alter the ordering of the cluster, we only need
## to match one of the columns.

expect_true(identical(scRep_fastF[, 1], scRep_fastT[, 1]) | identical(scRep_fastF[, 
    1], scRep_fastT[, 2]) | identical(scRep_fastF[, 1], scRep_fastT[, 3]))

expect_true(identical(scRep_fastF[, 2], scRep_fastT[, 1]) | identical(scRep_fastF[, 
    2], scRep_fastT[, 2]) | identical(scRep_fastF[, 2], scRep_fastT[, 3]))

expect_true(identical(scRep_fastF[, 3], scRep_fastT[, 1]) | identical(scRep_fastF[, 
    3], scRep_fastT[, 2]) | identical(scRep_fastF[, 3], scRep_fastT[, 3]))

#################### Testing if the parallel option yields identical result
#################### #################### mESC data is not large enough such that the
#################### parallelisation provides an improvement in time.  Ask YXL to test this
#################### code on an alterntive data.  t1 <- Sys.time() scRep_pT <- scReplicate( sce
#################### = sce_mESC, batch = sce_mESC$batch, kmeansK = c(1, 3, 3, 1, 1), fast_svd =
#################### TRUE, parallelParam = SnowParam(workers = 5, type = 'FORK')) t2 <-
#################### Sys.time() t3 <- Sys.time() scRep_pF <- scReplicate( sce = sce_mESC, batch
#################### = sce_mESC$batch, kmeansK = c(1, 3, 3, 1, 1), fast_svd = TRUE,
#################### parallelParam = BiocParallel::SerialParam()) t4 <- Sys.time()
#################### expect_identical(scRep_pF, scRep_pT) as.numeric(t2 - t1, units = 'secs')
#################### as.numeric(t4 - t3, units = 'secs') # expect_gt(as.numeric(t2 - t1, units
#################### = 'secs'), # as.numeric(t3 - t2, units = 'secs'))
