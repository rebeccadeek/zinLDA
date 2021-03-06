
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `zinLDA`: Zero-Inflated Latent Dirichlet Allocation

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/rebeccadeek/zinLDA.svg?branch=master)](https://travis-ci.com/rebeccadeek/zinLDA)
[![Codecov test
coverage](https://codecov.io/gh/rebeccadeek/zinLDA/branch/master/graph/badge.svg)](https://codecov.io/gh/rebeccadeek/zinLDA?branch=master)
<!-- badges: end -->

## About

Zero-inflated latent Dirichlet allocation is an unsupervised,
hierarchical, generative probabilistic model that facilitates
dimensionality reduction and detection of sparse latent clusters. This
package provides implementation of a Markov chain Monte Carlo (MCMC)
sampling procedure for the zinLDA model. Additionally, it provides a
method for simulating sparse count data from an underlying zinLDA model.

While the original paper developed this model for applications to
microbiome data and microbial subcommunity detection, it is flexible
enough to be used with numerous types of discrete count data.

## Installation

You can install the latest version of `zinLDA` from GitHub with:

``` r
install.packages("devtools")
devtools::install_github("rebeccadeek/zinLDA")
```

## Documentation and Examples

Help documentation for the `zinLDA` package is available in R. After
installing the package from GitHub via `devtools` and loading it with
`library()` use `?` to access the documentation for any of the four main
functions in the package. E.g.

``` r
?zinLDA
```

Additionally, the `zinLDA` package contains a vignette with a more
detailed description of the model, how it differs from currently
existing methods, and examples on how to simulate data from zinLDA and
fit the model. To include the vignette during installation and to access
it please use:

``` r
devtools::install_github("rebeccadeek/zinLDA", build_vignettes = TRUE)
vignette(package = "zinLDA")
```

Alternatively, there is a companion
[website](https://rebeccadeek.github.io/zinLDA/) that contains the same
statistical details and examples as the vignette.

## Contact

To report any bugs, issues, or suggestions please use the issue feature
on GitHub or contact the maintainer Rebecca Deek via
[email](mailto:rebecca.deek@pennmedicine.upenn.edu).
