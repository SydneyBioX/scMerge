# scMerge

`scMerge` is a R package for merging and normalising single-cell RNA-Seq datasets.


## Installation 

The installation process could take up to 5 minutes, depending if you have some of the packages pre-installed. 

``` r
# Some CRAN packages required by scMerge
install.packages(c("ruv", "rsvd", "igraph", "pdist", "proxy", "foreach", "doSNOW", "distr"))
devtools::install_github("theislab/kBET")

# Some BioConductor packages required by scMerge
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("SingleCellExperiment", "M3Drop"))


# Installing scMerge using
devtools::install_github("SydneyBioX/scMerge")
```

## Vignette

You can find the vignette at our website: https://sydneybiox.github.io/scMerge/index.html. 


## Reference

**scMerge: Integration of multiple single-cell transcriptomics datasets leveraging stable expression and pseudo-replication**

Yingxin Lin, Shila Ghazanfar, Kevin Y.X. Wang, Johann A. Gagnon-Bartsch, Kitty K. Lo, Xianbin Su, Ze-Guang Han, John T. Ormerod, Terence P. Speed, Pengyi Yang, Jean Y. H. Yang

doi: https://doi.org/10.1101/393280

BioRxiv preprint: https://www.biorxiv.org/content/early/2018/08/16/393280
