## ----style, echo = FALSE, results = 'asis'---------------------------------
BiocStyle::markdown()

## ----loading packages, warning = FALSE, message = FALSE--------------------
  suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scMerge)
  library(scater)
  })

## ----subsampling scMergeData, eval = FALSE, echo = FALSE-------------------
#  library(genefilter)
#  
#  load("~/Downloads/sce_mESC.rda")
#  data("segList_ensemblGeneID", package = "scMerge")
#  
#  set.seed(2019)
#  
#  example_sce = sce_mESC[, sce_mESC$batch %in% c("batch2", "batch3")]
#  example_sce$batch = droplevels(example_sce$batch)
#  batch2Sampled = sample(colnames(example_sce[,example_sce$batch == "batch2"]), 150)
#  batch3Sampled = sample(colnames(example_sce[,example_sce$batch == "batch3"]), 150)
#  
#  countsMat = SingleCellExperiment::counts(example_sce)
#  
#  batchTest = rowFtests(countsMat, fac = example_sce$batch)
#  celltypeTest = rowFtests(countsMat, fac = factor(example_sce$cellTypes))
#  
#  commonSegGenes = intersect(segList_ensemblGeneID$mouse$mouse_scSEG, rownames(sce_mESC))
#  
#  keepGenes = unique(c(commonSegGenes,
#                       rownames(batchTest)[rank(batchTest$p.value) < 200],
#                       rownames(celltypeTest)[rank(celltypeTest$p.value) < 1200]
#  ))
#  
#  example_sce = example_sce[keepGenes, c(batch2Sampled, batch3Sampled)]
#  example_sce = example_sce[rowSums(counts(example_sce)) != 0, colSums(counts(example_sce)) != 0]
#  
#  table(example_sce$batch,
#        example_sce$cellTypes)
#  
#  dim(example_sce)
#  
#  save(example_sce,
#       file = "data/example_sce.rda")

## ----loading data----------------------------------------------------------
## Subsetted mouse ESC data
data("example_sce", package = "scMerge")

## single-cell stably expressed gene list
data("segList_ensemblGeneID", package = "scMerge")

## ----checking raw data-----------------------------------------------------
scater::plotPCA(example_sce, 
                colour_by = "cellTypes", 
                shape_by = "batch",
                add_ticks = FALSE)

## ----unsupervised_default, results='hide',fig.show='hide'------------------
scMerge_unsupervised <- scMerge(sce_combine = example_sce, 
                                ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                kmeansK = c(3, 3),
                                assay_name = "scMerge_unsupervised")

## ----unsupervised_default_plotting-----------------------------------------
scater::plotPCA(scMerge_unsupervised, 
                colour_by = "cellTypes", 
                shape_by = "batch",
                run_args = list(exprs_values = "scMerge_unsupervised"),
                add_ticks = FALSE)

## ----unsupervised_prop1, results='hide',fig.show='hide'--------------------
scMerge_unsupervised_all <- scMerge(example_sce, 
                                    ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
                                    kmeansK = c(3, 3),
                                    assay_name = "scMerge_unsupervised_all",
                                    replicate_prop = 1)

## ----unsupervised_prop1_plotting-------------------------------------------
scater::plotPCA(scMerge_unsupervised_all, 
                colour_by = "cellTypes", 
                shape_by = "batch",
                run_args = list(exprs_values = "scMerge_unsupervised_all"),
                add_ticks = FALSE)

