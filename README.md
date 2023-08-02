
# scMerge

<br />

<img src="https://github.com/SydneyBioX/scMerge/raw/master/inst/logo.png" align="right" width="200" />

[![R build
status](https://github.com/SydneyBioX/scMerge/workflows/R-CMD-check/badge.svg)](https://github.com/SydneyBioX/scMerge/actions)
[![Codecov test
coverage](https://codecov.io/gh/SydneyBioX/scMerge/branch/master/graph/badge.svg)](https://codecov.io/gh/SydneyBioX/scMerge?branch=master)
[![](https://img.shields.io/badge/doi-10.1073/pnas.1820006116-blue.svg)](https://doi.org/10.1073/pnas.1820006116)
[![](https://img.shields.io/badge/devel%20version-1.5.0-blue.svg)](https://github.com/SydneyBioX/scMerge)
[![](https://img.shields.io/badge/download-1155/total-green.svg)](https://bioconductor.org/packages/stats/bioc/scMerge)
[![](https://img.shields.io/github/last-commit/SydneyBioX/scMerge.svg)](https://github.com/SydneyBioX/scMerge/commits/master)
[![](https://img.shields.io/badge/Docker%20image-available-blue.svg)](https://hub.docker.com/repository/docker/kevinwang09/scmerge)

<br />

`scMerge` is a R package for merging and normalising single-cell RNA-Seq
datasets.

## Installation

`scMerge` is available on Bioconductor
(<https://bioconductor.org/packages/scMerge>). You can install it using:

``` r
## Install scMerge from Bioconductor, requires R 3.6.0 or above
BiocManager::install("scMerge")
## You can also try to install the Bioconductor devel version of scMerge:
BiocManager::install("scMerge", version = "devel")
```

## Vignette

You can find the vignette at our website:

1. scMerge: <https://sydneybiox.github.io/scMerge/articles/scMerge.html>.
2. scMerge2: <https://sydneybiox.github.io/scMerge/articles/scMerge2.html>.

## Stably Expressed Genes

Stably expressed genes identified from this study can be extracted by

```
library(scMerge)
data(segList)
segList$human$human_scSEG # human SEG
segList$mouse$mouse_scSEG # mouse SEG
```

Or download csv files here (human SEG: [link](https://www.maths.usyd.edu.au/u/yingxinl/wwwnb/SEG/human_scSEG.csv); mouse SEG: [link](https://www.maths.usyd.edu.au/u/yingxinl/wwwnb/SEG/mouse_scSEG.csv))

For more detailed information and evaluation about SEG, please see our publication https://doi.org/10.1093/gigascience/giz106.

## Contact us

If you have any enquiries, especially about performing `scMerge`
integration on your own data, then please contact
<yingxin.lin@sydney.edu.au>. You can also [open an
issue](https://github.com/SydneyBioX/scMerge/issues) on GitHub.

## Reference


1. scMerge: **scMerge leverages factor analysis, stable expression, and
pseudoreplication to merge multiple single-cell RNA-seq datasets**. Yingxin Lin, Shila Ghazanfar, Kevin Y.X. Wang, Johann A. Gagnon-Bartsch,
Kitty K. Lo, Xianbin Su, Ze-Guang Han, John T. Ormerod, Terence P.
Speed, Pengyi Yang, Jean Y. H. Yang. (2019). Our manuscript published at PNAS can be found
[here](http://www.pnas.org/lookup/doi/10.1073/pnas.1820006116).

2. scMerge2: **Atlas-scale single-cell multi-sample multi-condition data integration using scMerge2**. Yingxin Lin, Yue Cao, Elijah Willie, Ellis Patrick, Jean Y.H. Yang. (2023). Our manuscript published in Nature Communications can be found [here](https://doi.org/10.1038/s41467-023-39923-2).
