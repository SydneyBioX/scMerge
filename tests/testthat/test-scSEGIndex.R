context("Test scSEGIndex")
data('example_sce', package = 'scMerge')
exprsMat = SummarizedExperiment::assay(example_sce, 'counts')[1:110, 1:20]
set.seed(1)
result = scSEGIndex(exprsMat = exprsMat, 
           cell_type = example_sce$cellTypes[1:20], ncore = 1)
gammaNormMix(data = exprsMat[1,])
