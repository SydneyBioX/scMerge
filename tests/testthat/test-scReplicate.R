context("Test scReplicate")
data("example_sce", package = "scMerge")

set.seed(100)

## Error case:
expect_error(
    scReplicate(sce_combine = example_sce, 
                batch = example_sce$batch, 
                replicate_prop = 1, 
                kmeansK = c(3, 3), 
                fast_svd = FALSE, 
                cell_type = example_sce$cellTypes, 
                cell_type_inc = 1:10, cell_type_match = TRUE))


## Case 1: supervised replicates: 
## Supply cell_type 
## Do not supply cell_type_inc 
## Do not upply cell_type_match
case1 = scReplicate(
    sce_combine = example_sce, 
    batch = example_sce$batch, 
    replicate_prop = 1, 
    kmeansK = c(3, 3), 
    fast_svd = FALSE, 
    cell_type = example_sce$cellTypes, 
    cell_type_inc = NULL, 
    cell_type_match = FALSE)

## Case 2.1: semi-supervised replicates version 2: 
## Supply cell_type 
## Do not supply cell_type_inc 
## Perform cell_type_match
## Do not supply marker 
## Do not supply marker_list
case2.1 = scReplicate(
    sce_combine = example_sce,
    batch = example_sce$batch,
    replicate_prop = 1, 
    kmeansK = c(3, 3), 
    fast_svd = FALSE, 
    cell_type = example_sce$cellTypes, 
    cell_type_inc = NULL, 
    cell_type_match = TRUE, 
    marker = NULL, 
    marker_list = NULL)


## Case 2.2: semi-supervised replicates version 2: 
## Supply cell_type 
## Do not supply cell_type_inc 
## Do not upply cell_type_match 
## Do not supply marker 
## Supply marker_list
case2.2 = scReplicate(
    sce_combine = example_sce, 
    batch = example_sce$batch,
    replicate_prop = 1,
    kmeansK = c(3, 3),
    fast_svd = FALSE, 
    cell_type = example_sce$cellTypes, 
    cell_type_inc = NULL, 
    cell_type_match = TRUE, 
    marker = NULL, 
    marker_list = 
        list(rownames(example_sce[1:10, ]), rownames(example_sce[1:10, ]))
)

## Case 2.3.1: semi-supervised replicates version 2: 
## Supply cell_type 
## Do not supply cell_type_inc 
## Do not supply cell_type_match
## Supply marker 
## Do not supply marker_list
case2.3.1 = scReplicate(sce_combine = example_sce, batch = example_sce$batch, replicate_prop = 1, kmeansK = c(3, 3), fast_svd = FALSE, cell_type = example_sce$cellTypes, 
                        cell_type_inc = NULL, cell_type_match = TRUE, marker = rownames(example_sce[1:10, ]), marker_list = NULL)

## Case 2.3.2: semi-supervised replicates version 2: 
## Supply cell_type 
## Do not supply cell_type_inc 
## Do not upply cell_type_match 
## Supply marker 
## Supply marker_list
case2.3.2 = scReplicate(sce_combine = example_sce, batch = example_sce$batch, replicate_prop = 1, kmeansK = c(3, 3), fast_svd = FALSE, cell_type = example_sce$cellTypes, 
                        cell_type_inc = NULL, cell_type_match = TRUE, marker = rownames(example_sce[1:10, ]), marker_list = list(rownames(example_sce[1:10, ]), 
                                                                                                                                 rownames(example_sce[1:10, ])))

## Because case2.3.1 and 2.3.2 do not depend on the input of marker_list
## so their output should be identical
expect_identical(case2.3.1, case2.3.2)



## Semi-supervised version 1 in the vignette
case3 = scReplicate(sce_combine = example_sce, batch = example_sce$batch, replicate_prop = 1, kmeansK = c(3, 3), fast_svd = FALSE, cell_type = example_sce$cellTypes, 
                    cell_type_inc = 1:10, cell_type_match = FALSE)

## Case 4.1: Do not Supply cell_type 
## Do not supply cell_type_inc 
## Do not supply cell_type_match 
## Do not supply marker 
## Do not supply marker_list
case4.1 = scReplicate(sce_combine = example_sce, batch = example_sce$batch, replicate_prop = 1, kmeansK = c(3, 3), fast_svd = FALSE, cell_type = NULL, 
                      cell_type_inc = NULL, cell_type_match = FALSE, marker = NULL, marker_list = NULL)


## Case 4.2: 
## Do not supply cell_type 
## Do not supply cell_type_inc 
## Do not supply cell_type_match 
## Do not supply marker 
## Supply marker_list
case4.2 = scReplicate(sce_combine = example_sce, batch = example_sce$batch, replicate_prop = 1, kmeansK = c(3, 3), fast_svd = FALSE, cell_type = NULL, 
                      cell_type_inc = NULL, cell_type_match = FALSE, marker = NULL, marker_list = list(rownames(example_sce[1:10, ]), rownames(example_sce[1:10, 
                                                                                                                                                           ])))

## Case 4.3: 
## Do not Supply cell_type 
## Do not supply cell_type_inc 
## Do not supply cell_type_match 
## Supply marker 
## Do not supply marker_list
case4.3 = scReplicate(sce_combine = example_sce, batch = example_sce$batch, replicate_prop = 1, kmeansK = c(3, 3), fast_svd = FALSE, cell_type = NULL, 
                      cell_type_inc = NULL, cell_type_match = FALSE, marker = rownames(example_sce[1:10, ]), marker_list = NULL)
