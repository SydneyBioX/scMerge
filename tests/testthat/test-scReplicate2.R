# More cases to test scReplicate
testthat::context("More cases to test scReplicate")

library(SingleCellExperiment)
library(scater)
library(scMerge)

# Test replicate output 

##################################################
##### Case that only has one cell type match #####
##################################################

subset_idx <- which((example_sce$cellTypes %in% c("2i", "a2i", "lif") & example_sce$batch %in% c("batch2")) |
                 (example_sce$cellTypes %in% c("a2i") & example_sce$batch %in% c("batch3")))


example_sce1 = example_sce[, subset_idx]

set.seed(2)
index = sample(seq_len(ncol(example_sce1)))

example_sce2 = example_sce1[, index]

repMat1 = scReplicate(
  sce_combine = example_sce1, 
  batch = example_sce1$batch, 
  replicate_prop = 1, 
  kmeansK = c(3, 1), 
  fast_svd = FALSE, 
  cell_type = NULL, 
  cell_type_inc = NULL, 
  cell_type_match = FALSE, 
  marker = NULL, 
  marker_list = NULL,
  verbose = TRUE,
  return_all = TRUE)

expect_equal(max(repMat1$mnc$group), 3)

repMat2 = scReplicate(
  sce_combine = example_sce2, 
  batch = example_sce2$batch, 
  replicate_prop = 1, 
  kmeansK = c(3, 1), 
  fast_svd = FALSE, 
  cell_type = NULL, 
  cell_type_inc = NULL, 
  cell_type_match = FALSE, 
  marker = NULL, 
  marker_list = NULL,
  verbose = TRUE,
  return_all = TRUE)


repMat2$repMat <- repMat2$repMat[order(index),]
repMat2$replicate_vector <- repMat2$replicate_vector[order(index)]
tab <- table(repMat2$replicate_vector, repMat1$replicate_vector)

expect_equal(sum(rowSums(tab == 0) == 2), 3)

##################################################
##### Case that have two cell type match     #####
##################################################


subset_idx <- which((example_sce$cellTypes %in% c("2i", "a2i", "lif") & example_sce$batch %in% c("batch2")) |
                      (example_sce$cellTypes %in% c("2i", "a2i") & example_sce$batch %in% c("batch3")))


example_sce1 = example_sce[, subset_idx]

set.seed(2)
index = sample(seq_len(ncol(example_sce1)))

example_sce2 = example_sce1[, index]

repMat1 = scReplicate(
  sce_combine = example_sce1, 
  batch = example_sce1$batch, 
  replicate_prop = 1, 
  kmeansK = c(3, 2), 
  fast_svd = FALSE, 
  cell_type = NULL, 
  cell_type_inc = NULL, 
  cell_type_match = FALSE, 
  marker = NULL, 
  marker_list = NULL,
  verbose = TRUE,
  return_all = TRUE)

expect_equal(max(repMat1$mnc$group), 3)

repMat2 = scReplicate(
  sce_combine = example_sce2, 
  batch = example_sce2$batch, 
  replicate_prop = 1, 
  kmeansK = c(3, 2), 
  fast_svd = FALSE, 
  cell_type = NULL, 
  cell_type_inc = NULL, 
  cell_type_match = FALSE, 
  marker = NULL, 
  marker_list = NULL,
  verbose = TRUE,
  return_all = TRUE)


repMat2$repMat <- repMat2$repMat[order(index),]
repMat2$replicate_vector <- repMat2$replicate_vector[order(index)]
tab <- table(repMat2$replicate_vector, repMat1$replicate_vector)

expect_equal(sum(rowSums(tab == 0) == 2), 3)


##################################################
#####       Case that have all ones          #####
##################################################

subset_idx <- which((example_sce$cellTypes %in% c("2i", "a2i") & example_sce$batch %in% c("batch2")) |
                      (example_sce$cellTypes %in% c("a2i") & example_sce$batch %in% c("batch3")))


example_sce1 = example_sce[, subset_idx]

example_sce1$batch <- as.character(example_sce1$batch)
example_sce1$batch[example_sce1$batch == "batch2" & example_sce1$cellTypes == "a2i"] <- "batch1"

set.seed(2)
index = sample(seq_len(ncol(example_sce1)))

example_sce2 = example_sce1[, index]

repMat1 = scReplicate(
  sce_combine = example_sce1, 
  batch = example_sce1$batch, 
  replicate_prop = 1, 
  kmeansK = c(1, 1, 1), 
  fast_svd = FALSE, 
  cell_type = NULL, 
  cell_type_inc = NULL, 
  cell_type_match = FALSE, 
  marker = NULL, 
  marker_list = NULL,
  verbose = TRUE,
  return_all = TRUE)

mnc_expect <- data.frame(group = c(1, 2, 1),
                         batch = c(1, 2, 3),
                         cluster = c(1, 1, 1))

expect_identical(repMat1$mnc, mnc_expect)

repMat2 = scReplicate(
  sce_combine = example_sce2, 
  batch = example_sce2$batch, 
  replicate_prop = 1, 
  kmeansK = c(1, 1, 1), 
  fast_svd = FALSE, 
  cell_type = NULL, 
  cell_type_inc = NULL, 
  cell_type_match = FALSE, 
  marker = NULL, 
  marker_list = NULL,
  verbose = TRUE,
  return_all = TRUE)


repMat2$repMat <- repMat2$repMat[order(index),]
repMat2$replicate_vector <- repMat2$replicate_vector[order(index)]

expect_identical(repMat2$mnc, mnc_expect)

expect_equivalent(repMat1$repMat, repMat2$repMat)

