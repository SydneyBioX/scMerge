# scMerge

`scMerge` is a R package for merging and normalising single-cell RNA-Seq datasets.

## Installation 

The installation process could take up to 5 minutes, depending if you have some of the packages pre-installed. 

``` r
# Some CRAN packages required by scMerge
install.packages(c("ruv", "rsvd", "igraph", "pdist", "proxy"))
devtools::install_github("theislab/kBET")

# Some BioConductor packages required by scMerge
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("SingleCellExperiment", "M3Drop"))


# Installing scMerge using
devtools::install_github("SydneyBioX/scMerge")
```




## Reference
Our paper is available on bioRxiv. 
