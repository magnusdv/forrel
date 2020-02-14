<!-- README.md is generated from README.Rmd. Please edit that file -->

forrel <img src="man/figures/logo.png" align="right" height=140/>
=================================================================

Introduction
------------

The goal of `forrel` is to provide forensic pedigree computations and
relatedness inference from genetic marker data. Some of its
functionality is derived from `paramlink`, which is no longer actively
developed. ‘forrel’ is part of the *ped suite*, a collection of packages
for pedigree analysis with `pedtools` as the central package for
handling pedigrees and marker data.

The most important analyses currently supported by forrel are:

-   Likelihood ratio (LR) computations for relationship testing
-   Pairwise relatedness inference: Estimation of IBD (kappa)
    coefficients from marker data
-   Visualisation of IBD coefficients in the IBD triangle
-   Simulation of marker genotypes. Unconditional or conditional on
    known genotypes
-   Power analysis for relationship testing: LR distributions, exclusion
    power (EP) and inclusion power (IP)
-   Tailor-made functions for power analysis in family reunion cases:
    -   `missingPersonPlot()`
    -   `missingPersonEP()`
    -   `missingPersonIP()`
    -   `MPPsims()`
    -   `powerPlot()`
-   Import of pedigree data and frequency databases from the `Familias`
    software.

Installation
------------

To get the current official version of forrel, install from CRAN as
follows:

``` r
install.packages("forrel")
```

Alternatively, you can obtain the latest development version from
GitHub:

``` r
# install.packages("devtools") # install devtools if needed
devtools::install_github("magnusdv/forrel")
```

An example
----------

In this short introduction, we first demonstrate simulation of marker
data for a pair of siblings. Then - pretending the relationship is
unknown to us - we estimate the relatedness between the brothers using
the simulated data. If all goes well, the estimate should be close to
the expected value for siblings.

``` r
library(forrel)
#> Loading required package: pedtools
```

**Create the pedigree**

We start by creating and plotting a pedigree with two brothers, named
`bro1` and `bro2`.

``` r
bros = c("bro1", "bro2")
x = nuclearPed(children = bros)
plot(x)
```

<img src="man/figures/README-sibs-1.png" style="display: block; margin: auto;" />

**Marker simulation**

Now let us simulate the genotypes of 100 independent SNPs for the
brothers. Each SNP has alleles 1 and 2, with equal frequencies by
default. This is an example of *unconditional* simulation, since we
don’t give any genotypes to condition on. Unconditional simulation is
performed by simple *gene dropping*, i.e., by drawing random alleles
independently for the parents, followed by a “Mendelian coin toss” in
each parent-child transmission.

``` r
x = markerSim(x, N = 100, ids = bros, alleles = 1:2, seed = 1234)
#> Unconditional simulation of 100 autosomal markers.
#> Individuals: bro1, bro2
#> Allele frequencies:
#>    1   2
#>  0.5 0.5
#> Mutation model: No 
#> 
#> Simulation finished.
#> Number of calls to the likelihood function: 0.
#> Total time used: 0.06 seconds.
```

Note 1: The `seed` argument is passed onto the random number generator.
If you use the same seed, you should get exactly the same results.  
Note 2: To suppress the informative messages printed during simulation,
add `verbose = FALSE` to the function call.

The pedigree `x` now has 100 markers attached to it. The genotypes of
the first few markers are shown when printing `x` to the screen:

``` r
x
#>    id fid mid sex <1> <2> <3> <4> <5>
#>     1   *   *   1 -/- -/- -/- -/- -/-
#>     2   *   *   2 -/- -/- -/- -/- -/-
#>  bro1   1   2   1 1/1 1/2 1/1 2/1 2/2
#>  bro2   1   2   1 1/1 2/1 1/1 2/1 2/2
#> Only 5 (out of 100) markers are shown. See `?print.ped` for options.
```

**Estimation of IBD coefficients**

The `IBDestimate()` function estimates the coefficients of
*identity-by-descent* (IBD) between pairs of individuals, from the
available marker data.

``` r
k = IBDestimate(x, ids = bros)
k
#>    ID1  ID2   N        k0        k1        k2
#> 1 bro1 bro2 100 0.1484924 0.5513491 0.3001585
```

The theoretical expectation for non-inbred full siblings is
(*κ*<sub>0</sub>, *κ*<sub>1</sub>, *κ*<sub>2</sub>) = (0.25, 0.5, 0.25).
To get a visual sense of how close our estimate is, it is instructive to
plot it in the IBD triangle:

``` r
showInTriangle(k, labels = TRUE)
```

<img src="man/figures/README-triangle-1.png" style="display: block; margin: auto;" />
