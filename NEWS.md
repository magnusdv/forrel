# forrel 1.3.0

# New features

* The new function `ibdEstimate()` replaces the previous `IBDestimate()` (note the name change). This is a complete rewrite, which optimises the log-likelihood using a projected gradient descent algorithm, combined with a version of Armijo line search.

* The function `ibdBootstrap()` replaces the previous `kappaBootstrap()` and `deltaBootstrap()`, and is considerably faster. This function implements both parametric and non-parametric bootstrap, controlled with the `method` parameter.

* The output of `ibdEstimate()` now has a class attribute "ibdEst", with its own print and subsetting methods.



# forrel 1.2.0

## New features

* `kinshipLR()` now handles linked markers by wrapping MERLIN.

* New functions `kappaBootstrap()` and `deltaBootstrap()` for assessing the uncertainty of pairwise relatedness estimates.

* New function `randomPersonEP()` handling a common special case of `exclusionPower()`.

## Other changes

* forrel now depends on version 0.9.6 (or later) of pedtools.

* Deprecated arguments `id.labels` and `frametitles` in `missingPersonPlot()` has been removed.


# forrel 1.1.0

## New features

* Implement parallelisation in `profileSim()`.

* Partial rewrite of `kinshipLR()`, including new argument `source`.

* Added the `NorwegianFrequencies` dataset, containing allele frequencies for 35 STR markers.

* New function `missingPersonLR()`.

* New function `checkPairwise()` replaces the (long obsolete) `examineKinships()`.

* New functions `markerSimParametric()` and `profileSimParametric()` for simulating marker data for two individuals with given kappa (or condensed identity) coefficients.

## Bug fixes

* In `profileSim()`, fix bug resulting in identical seeds given to each parallel cluster.


# forrel 1.0.1

## New features

* `readFam()` now has a parameter `Xchrom` which can be used to indicate that the markers included in the file are on the X chromosome

* `MPPsims()` is more flexible, and allows subsetting of its output.

* `powerPlot()` is more flexible and allows finer control of the plot contents

## Bug fixes

* Fixed several glitches in `readFam()`. It is more robust now, and fails gracefully in certain situations which cannot currently be handled (e.g. if the file contains twins).


# forrel 1.0.0

* Initial CRAN release
