[![Travis build status](https://travis-ci.org/SydneyBioX/scMerge.svg?branch=master)](https://travis-ci.org/SydneyBioX/scMerge)

# scMerge

`scMerge` is a R package for merging and normalising single-cell RNA-Seq datasets.


## Installation 

The installation process could take up to 5 minutes, depending if you have some of the packages pre-installed. 

``` r
## Some CRAN packages required by scMerge
install.packages(c("BiocManager", "cluster", "distr", "doSNOW", "foreach", "igraph", "irlba", "iterators", "matrixStats", "pdist", "proxy",  "Rcpp", "RcppEigen", "rsvd", "ruv"))

## Some BioConductor packages required by scMerge
BiocManager::install(c("BiocParallel", "M3Drop", "SingleCellExperiment"))

## Installing scMerge on R 3.6
devtools::install_github("SydneyBioX/scMerge")

## Installing scMerge on R 3.5
devtools::install_github("SydneyBioX/scMerge", ref = "280a7724cd49f10e535f4023a085bc6dc7dc432d")
```

## Vignette

You can find the vignette at our website: https://sydneybiox.github.io/scMerge/index.html. 


## Case studies

You can find a list of case studies here: https://sydneybiox.github.io/scMerge/articles/.


## Contact us

If you have any enquires, especially about performing `scMerge` integration on your data, then please contact bioinformatics@maths.usyd.edu.au. 

## Reference

**scMerge: Integration of multiple single-cell transcriptomics datasets leveraging stable expression and pseudo-replication**

Yingxin Lin, Shila Ghazanfar, Kevin Y.X. Wang, Johann A. Gagnon-Bartsch, Kitty K. Lo, Xianbin Su, Ze-Guang Han, John T. Ormerod, Terence P. Speed, Pengyi Yang, Jean Y. H. Yang

doi: https://doi.org/10.1101/393280

BioRxiv preprint: https://www.biorxiv.org/content/early/2018/08/16/393280
