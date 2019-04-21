
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMerge <img src="./inst/logo.png" align="right" width="200" />

[![Travis build
status](https://travis-ci.org/SydneyBioX/scMerge.svg?branch=master)](https://travis-ci.org/SydneyBioX/scMerge)
[![Codecov test
coverage](https://codecov.io/gh/SydneyBioX/scMerge/branch/master/graph/badge.svg)](https://codecov.io/gh/SydneyBioX/scMerge?branch=master)

`scMerge` is a R package for merging and normalising single-cell RNA-Seq
datasets.

## Installation

The installation process could take up to 5 minutes, depending if you
have some of the packages pre-installed.

``` r
## Some CRAN packages required by scMerge
install.packages(c("BiocManager", "cluster", "distr", "doSNOW", 
                   "foreach", "igraph", "irlba", "iterators", 
                   "matrixStats", "pdist", "proxy",  "Rcpp", 
                   "RcppEigen", "rsvd", "ruv"))

## Some BioConductor packages required by scMerge
BiocManager::install(c("BiocParallel", "M3Drop", "SingleCellExperiment"))

## Install scMerge from BioConductor, requires R 3.6.0 or above
BiocManager::install("scMerge", version = "3.9")


## Install the development version of scMerge
install_github("SydneyBioX/scMerge", 
               build_opts = c("--no-resave-data", "--no-manual"),
               build_vignettes = TRUE)
```

## Vignette

You can find the vignette at our website:
<https://sydneybiox.github.io/scMerge/index.html>.

## Case studies

You can find a list of case studies here:
<https://sydneybiox.github.io/scMerge/articles/>.

## Contact us

If you have any enquries, especially about performing `scMerge`
integration on your data, then please contact
<bioinformatics@maths.usyd.edu.au>.

## Reference

**scMerge: Integration of multiple single-cell transcriptomics datasets
leveraging stable expression and pseudo-replication**

Yingxin Lin, Shila Ghazanfar, Kevin Y.X. Wang, Johann A. Gagnon-Bartsch,
Kitty K. Lo, Xianbin Su, Ze-Guang Han, John T. Ormerod, Terence P.
Speed, Pengyi Yang, Jean Y. H. Yang

doi: <https://doi.org/10.1101/393280>

BioRxiv preprint:
<https://www.biorxiv.org/content/early/2018/08/16/393280>
