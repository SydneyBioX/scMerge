
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMerge

<br />

<img src="https://github.com/SydneyBioX/scMerge/raw/master/inst/logo.png" align="right" width="200" />

[![Travis build
status](https://travis-ci.org/SydneyBioX/scMerge.svg?branch=master)](https://travis-ci.org/SydneyBioX/scMerge)
[![Codecov test
coverage](https://codecov.io/gh/SydneyBioX/scMerge/branch/master/graph/badge.svg)](https://codecov.io/gh/SydneyBioX/scMerge?branch=master)
[![](https://img.shields.io/badge/doi-10.1073/pnas.1820006116-blue.svg)](https://doi.org/10.1073/pnas.1820006116)
[![](https://img.shields.io/badge/devel%20version-1.1.2-blue.svg)](https://github.com/SydneyBioX/scMerge)
[![](https://img.shields.io/badge/download-84/total-green.svg)](https://bioconductor.org/packages/stats/bioc/scMerge)
[![](https://img.shields.io/github/last-commit/SydneyBioX/scMerge.svg)](https://github.com/SydneyBioX/scMerge/commits/master)

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
<https://sydneybiox.github.io/scMerge/index.html>.

## Case studies

You can find a list of case studies here:
<https://sydneybiox.github.io/scMerge/articles/>.

## Contact us

If you have any enquries, especially about performing `scMerge`
integration on your own data, then please contact
<bioinformatics@maths.usyd.edu.au>. You can also [open an
issue](https://github.com/SydneyBioX/scMerge/issues) on GitHub.

## Reference

**scMerge leverages factor analysis, stable expression, and
pseudoreplication to merge multiple single-cell RNA-seq datasets**

Yingxin Lin, Shila Ghazanfar, Kevin Y.X. Wang, Johann A. Gagnon-Bartsch,
Kitty K. Lo, Xianbin Su, Ze-Guang Han, John T. Ormerod, Terence P.
Speed, Pengyi Yang, Jean Y. H. Yang

Our manuscript published at PNAS can be found
[here](http://www.pnas.org/lookup/doi/10.1073/pnas.1820006116).
