# scMerge 1.1.6

* Accepts `DelayedArray`, `HDF5Array` and `dgCMatrix` inputs in the slots of input SCE objects. 
* Significant speed optimisation on `scSEGIndex` and add `BiocParallel` support. 
* Updated `scSEGIndex` references after publication. 
* `scMerge` now has the `svd_k` input that controls the number eigenvectors needed in the RUV step to allow fast approximation for large dataset. 
* Now using `BiocSingular` to manage all SVD components. 
* Now automatically remove zeroes in the rows and columns of the SCE. 

# scMerge 1.1.5

* Adding version restrictions on `S4Vectors` and `SingleCellExperiment` dependent packages. 

# scMerge 1.1.4

* `plot_igraph` would allow suppression of igraph output during unsupervised scMerge

# scMerge 1.1.3

* Column name must be non-NULL and without duplicates

# scMerge 1.1.2

* Resolved problems with only a single linking cell-type across multiple batches

# scMerge 1.1.0

* Accepted by Bioconductor

# scMerge 0.99 (development version)

## scMerge 0.99.24
* Updated ciation information due to PNAS acceptance. 

## scMerge 0.99.23
* Fixed assignment based on feedbacks

## scMerge 0.99.21
* Increase code coverage to 85%.


## scMerge 0.99.20
* Updated vignette on SEGs and manuals


## scMerge 0.99.19
* Fixed spelling
* Added verbose option
* Code coverage at 75 percent (more tests on error handling needed)


## scMerge 0.99.17 
* Fixed README `install_github` vignette issue. 
* Fixed pkgdown organisation issue.
* Major updates on the scReplicate function: more informative messages. 
* Using cross-product of matrix to perform SVD to speed up calculations.
* Added testing scripts. 
* Fixed vignette text output issues.

## scMerge 0.99.11 
* Reduced data size in scMerge to pass BioC checks