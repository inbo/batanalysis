
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![Project Status: Concept – Minimal or no implementation has been done
yet, or the repository is only intended to be a limited example, demo,
or
proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License](https://img.shields.io/badge/license-GPL--3-blue.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![Release](https://img.shields.io/github/release/inbo/batanalysis.svg)](https://github.com/inbo/batanalysis/releases)
![GitHub](https://img.shields.io/github/license/inbo/batanalysis) [![R
build
status](https://github.com/inbo/batanalysis/workflows/check%20package%20on%20main/badge.svg)](https://github.com/inbo/batanalysis/actions)
![r-universe
name](https://inbo.r-universe.dev/badges/:name?color=c04384)
![r-universe package](https://inbo.r-universe.dev/badges/batanalysis)
[![Codecov test
coverage](https://codecov.io/gh/inbo/batanalysis/branch/main/graph/badge.svg)](https://app.codecov.io/gh/inbo/batanalysis?branch=main)
![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/inbo/batanalysis.svg)
![GitHub repo
size](https://img.shields.io/github/repo-size/inbo/batanalysis.svg)
<!-- badges: end -->

# `batanalysis`

The goal of `batanalysis` is to analyse the monitoring data from bats in
Flanders.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("inbo/batanalysis")
```

## Example

``` r
library(batanalysis)
# connection to the database
constring <- "replace this with the correct connection string"
origin <- odbc::dbConnect(odbc::odbc(), .connection_string = constring)
# folder to store the imported data
target <- tempfile()
dir.create(target)
target <- git2r::init(target)
# import the data
import_raw_data(origin = origin, target = target)
```
