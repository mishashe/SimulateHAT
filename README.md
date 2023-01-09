SimulateHAT
================
2023-01-09

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package description

<!-- badges: start -->
<!-- badges: end -->

The goal of SimulateHAT is to segment genome with mutations to segments
with significantly different mutations density.

## Installation

You can install the development version of SimulateHAT from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("mishashe/SimulateHAT")
```

To load necessary libraries:

``` r
library(SimulateHAT)
library(RColorBrewer, quietly=T) # to plot results
library(stringr)
library(tidyverse)
library(scales)
library(ggExtra)
library(stringi)
library(ggplot2)
library(plotrix)
library(Rcpp)
library(RcppArmadillo)
library(roxygen2)
options(warn=-1)
```
