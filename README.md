<!-- README.md is generated from README.Rmd. Please edit that file -->
forrel <img src="man/figures/logo.png" align="right" height=140/>
=================================================================

Introduction
------------

The goal of `forrel` is to provide forensic pedigree computations and relatedness estimation from genetic marker data. Some of its functionality is derived from the retired `paramlink` package, but updated to comply with the `pedtools` package for handling pedigrees and marker data.

`forrel` is under development, and only contains a few functions for now. The most important of these are:

-   `markerSim()` : Simulation of marker data
-   `Familias2ped()` : Conversion of `Familias` objects to `pedtools::ped` objects

Installation
------------

To get the lastest version, install from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install forrel from GitHub
devtools::install_github("magnusdv/forrel")
```
