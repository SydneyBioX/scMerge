context("Test scMerge2h")


## Loading example data
data('example_sce', package = 'scMerge')
data('segList_ensemblGeneID', package = 'scMerge')



# Create sample information
example_sce$sample <- rep(c(1:4), each = 50)

example_sce$sample2 <- c(sample(c(1:2), 50, replace = TRUE),
                         sample(c(3:4), 50, replace = TRUE),
                         sample(c(5:6), 50, replace = TRUE),
                         sample(c(7:8), 50, replace = TRUE))


#############################################################################
##### Case 1: Level 1: sample level, Level 2: dataset level #################
#############################################################################

h_idx_list <- list(level1 = split(1:200, example_sce$batch),
                   level2 = list(1:200))



batch_list <- list(level1 = split(example_sce$sample, example_sce$batch),
                   level2 = list(example_sce$batch))





level1_subset_idx1 <- h_idx_list$level1[[1]]
level1_batch_idx1 <- batch_list$level1[[1]]
exprsMat_level1_idx1 <- scMerge2(exprsMat = logcounts(example_sce[, level1_subset_idx1]),
                                 batch = level1_batch_idx1,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2)
level1_subset_idx2 <- h_idx_list$level1[[2]]
level1_batch_idx2 <- batch_list$level1[[2]]
exprsMat_level1_idx2 <- scMerge2(exprsMat = logcounts(example_sce[, level1_subset_idx2]),
                                 batch = level1_batch_idx2,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2)


input_mat <- cbind(exprsMat_level1_idx1$newY,
                   exprsMat_level1_idx2$newY)
level2_batch <- batch_list$level2[[1]]
exprsMat_level2 <- scMerge2(exprsMat = input_mat,
                            batch = level2_batch,
                            ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                            ruvK = 5,
                            chosen.hvg = rownames(input_mat))


exprsMat <- scMerge2h(exprsMat = logcounts(example_sce),
                      batch_list = batch_list,
                      h_idx_list = h_idx_list,
                      ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                      ruvK_list = c(2, 5))




expect_equal(as(exprsMat[[1]][, colnames( exprsMat_level1_idx1$newY)], "dgCMatrix"),
             as(exprsMat_level1_idx1$newY, "dgCMatrix"))


expect_equal(as(exprsMat[[1]][, colnames( exprsMat_level1_idx2$newY)], "dgCMatrix"),
             as(exprsMat_level1_idx2$newY, "dgCMatrix"))



expect_equal(as(exprsMat[[2]][, colnames( exprsMat_level2$newY)], "dgCMatrix"),
             as(exprsMat_level2$newY, "dgCMatrix"))








#############################################################################
## Case 2: Level 1: sample level (only one batch), Level 2: dataset level ###
#############################################################################

h_idx_list <- list(level1 = split(1:200, example_sce$batch),
                   level2 = list(1:200))



batch_list <- list(level1 = split(example_sce$sample, example_sce$batch),
                   level2 = list(example_sce$batch))


h_idx_list$level1$batch3 <- NULL
batch_list$level1$batch3 <- NULL

level1_subset_idx1 <- h_idx_list$level1[[1]]
level1_batch_idx1 <- batch_list$level1[[1]]
exprsMat_level1_idx1 <- scMerge2(exprsMat = logcounts(example_sce[, level1_subset_idx1]),
                                 batch = level1_batch_idx1,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2)
level1_subset_idx2 <-setdiff(h_idx_list$level2[[1]], h_idx_list$level1[[1]])


input_mat <- cbind(exprsMat_level1_idx1$newY,
                   DelayedArray(logcounts(example_sce)[, level1_subset_idx2]))
level2_batch <- batch_list$level2[[1]]
exprsMat_level2 <- scMerge2(exprsMat = input_mat,
                            batch = level2_batch,
                            ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                            ruvK = 5,
                            chosen.hvg = rownames(input_mat))


exprsMat <- scMerge2h(exprsMat = logcounts(example_sce),
                      batch_list = batch_list,
                      h_idx_list = h_idx_list,
                      ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                      ruvK_list = c(2, 5))




expect_equal(as(exprsMat[[1]][, colnames( exprsMat_level1_idx1$newY)], "dgCMatrix"),
             as(exprsMat_level1_idx1$newY, "dgCMatrix"))



expect_equal(as(exprsMat[[2]][, colnames( exprsMat_level2$newY)], "dgCMatrix"),
             as(exprsMat_level2$newY, "dgCMatrix"))




###################################################################################################
## Case 3: Level 1: sample level, Level 2: dataset level (skip one batch), Level 3: study level ###
###################################################################################################

h_idx_list <- list(level1 = split(1:ncol(example_sce), example_sce$sample),
                   level2 = split(1:ncol(example_sce), example_sce$batch),
                   level3 = list(1:ncol(example_sce)))



batch_list <- list(level1 = split(example_sce$sample2, example_sce$sample),
                   level2 = split(example_sce$sample, example_sce$batch),
                   level3 = list(example_sce$batch))


h_idx_list$level2[[2]] <- NULL
batch_list$level2[[2]] <- NULL

level1_subset_idx1 <- h_idx_list$level1[[1]]
level1_batch_idx1 <- batch_list$level1[[1]]
exprsMat_level1_idx1 <- scMerge2(exprsMat = logcounts(example_sce[, level1_subset_idx1]),
                                 batch = level1_batch_idx1,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2)

level1_subset_idx2 <- h_idx_list$level1[[2]]
level1_batch_idx2 <- batch_list$level1[[2]]
exprsMat_level1_idx2 <- scMerge2(exprsMat = logcounts(example_sce[, level1_subset_idx2]),
                                 batch = level1_batch_idx2,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2)

level1_subset_idx3 <- h_idx_list$level1[[3]]
level1_batch_idx3 <- batch_list$level1[[3]]
exprsMat_level1_idx3 <- scMerge2(exprsMat = logcounts(example_sce[, level1_subset_idx3]),
                                 batch = level1_batch_idx3,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2)



level1_subset_idx4 <- h_idx_list$level1[[4]]
level1_batch_idx4 <- batch_list$level1[[4]]
exprsMat_level1_idx4 <- scMerge2(exprsMat = logcounts(example_sce[, level1_subset_idx4]),
                                 batch = level1_batch_idx4,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2)


input_mat <- cbind(exprsMat_level1_idx1$newY,
                   exprsMat_level1_idx2$newY)
level2_batch_idx1 <- batch_list$level2[[1]]
exprsMat_level2_idx1 <- scMerge2(exprsMat = input_mat,
                                 batch = level2_batch_idx1,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 2,
                                 chosen.hvg = rownames(input_mat))



input_mat <- cbind(exprsMat_level2_idx1$newY,
                   cbind(exprsMat_level1_idx3$newY,
                         exprsMat_level1_idx4$newY))
level3_batch_idx1 <- batch_list$level3[[1]]
input_mat <- input_mat[, colnames(example_sce)]
exprsMat_level3 <- scMerge2(exprsMat = input_mat,
                                 batch = level3_batch_idx1,
                                 ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                 ruvK = 5,
                                 chosen.hvg = rownames(input_mat))



exprsMat <- scMerge2h(exprsMat = logcounts(example_sce),
                      batch_list = batch_list,
                      h_idx_list = h_idx_list,
                      ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                      ruvK_list = c(2, 2, 5))




expect_equal(as(exprsMat[[1]][, colnames( exprsMat_level1_idx1$newY)], "dgCMatrix"),
             as(exprsMat_level1_idx1$newY, "dgCMatrix"))


expect_equal(as(exprsMat[[1]][, colnames( exprsMat_level1_idx2$newY)], "dgCMatrix"),
             as(exprsMat_level1_idx2$newY, "dgCMatrix"))


expect_equal(as(exprsMat[[1]][, colnames( exprsMat_level1_idx3$newY)], "dgCMatrix"),
             as(exprsMat_level1_idx3$newY, "dgCMatrix"))


expect_equal(as(exprsMat[[2]][, colnames( exprsMat_level2_idx1$newY)], "dgCMatrix"),
             as(exprsMat_level2_idx1$newY, "dgCMatrix"))




expect_equal(as(exprsMat[[3]][, colnames( exprsMat_level3$newY)], "dgCMatrix"),
             as(exprsMat_level3$newY, "dgCMatrix"))










