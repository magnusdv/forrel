# forrel 1.7.0

## New features

* New function `quickLR()` performs the most common kinship tests (paternity and sibship) for a pair of individuals. 

* New function `ibdLoglik()` computes the pairwise likelihood for given set IBD coefficients (kappa or delta). It replaces the previous (unexported) `.IBDlikelihood` and versions of this.

* `checkPairwise()` gains argument `plotType`, which can be either value is either "base" (default), "ggplot" or "plotly". 

* In `checkPairwise()` the output column LR is renamed to GLR (generalised likelihood ratio).

* `checkPairwise()` can now estimate p-values for each pairwise relationship, by simulating the distribution of GLR under the null. This feature is controlled through new arguments `nsim`, `pvalThreshold` and `seed`. By default `nsim = 0`, meaning that no simulations are performed.

* `checkPairwise()` gains argument `ids`, allowing to restrict checks to those individuals.

* `checkPairwise()` now contains a detailed description of each relationship according to the pedigree, obtained with the **verbalisr** package.

* Improved progress bars using the `pbapply` package.

## Other

* Refresh code in examples.

* Removed (long deprecated) `IBDtriangle()`, replaced by `ribd::ibdTriangle()`.

* Improved `missingPersonPlot()` appearance.

* Better removal of missing data in `ibdEstimate()`.


# forrel 1.6.1

* Functions facilitating interaction with the Familias software have been moved to a dedicated package, **pedFamilias**. However, `readFam()` and `writeFam()` will continue to be re-exported from **forrel** for some time.

* Remove dependency of **pedmut**.

* In `kinshipLR()`, improve checking of argument `markers`.


# forrel 1.6.0

## New features

* New function `writeFam()`.

* `readFam()` gains argument `fallbackModel` to be applied if the function encounters unrecognisable mutation models in the .fam file.

## Other

* Remove use of `packageVersion()` (thereby sidestepping CRAN issue).

* Use `pedtools::setMutmod()` instead of `pedprobr::setMutationModel()`.

* Update version requirements for imported packages.


# forrel 1.5.2

## New features

* `readFam()` includes more input checks and gives better error/warning messages.


# forrel 1.5.1

## New features

* `readFam()` now handle all mutation models available in the Familias software.

* `readFam()` and `Familias2ped()` gain argument `prefixAdded`.

## Bug fixes

* Fixed bug in `missingPersonLR()` manifesting when POI and MP have identical names.

## Other

* Added license statement GPL (>=2).

* Expanded README, adding an example with conditional simulation.

* Updated package documentation.

* Fixed CRAN note by avoiding `ibdsim2:::`.


# forrel 1.5.0

## Breaking change
* The default output of `profileSim(x, N = 1)` is now a single pedigree, instead of a list of length 1 (containing the pedigree). This is usually the desired output of `N = 1` in interactive use. To override this behaviour, set `simplify1 = FALSE`.

## New features

* Added dataset FORCE describing the FORCE snp panel (Tillmar et al, 2021, doi:10.3390/genes12121968).

* `profileSim()` has a new argument `simplify1` (by default TRUE) controlling the output when `N = 1`.

* `profileSim()` now allows `markers` to be a list of frequency vectors, simplifying the code in unconditional simulations. For instance, the following command now works: `nuclearPed() |> profileSim(markers = NorwegianFrequencies)`. Previously this required an intermediate step of `setMarkers()`.

* In `kinshipLR()` the treatment of linked markers (with MERLIN) has been rewritten and is substantially more efficient. A new argument `keepMerlin` allows to keep merlin files for debugging.

* `missingPersonLR()` was overhauled, making it more user friendly.

* `missingPersonPlot()` has been modified and updated in sync with changes in **pedtools**.


# forrel 1.4.1

* Fix CRAN complaint: Use `inherits()` instead of `class()`.


# forrel 1.4.0

## New features

* New function `findExclusions()` for identifying incompatible markers in identification cases.

* `powerPlot()` gains a logical argument `jitter`, which can be switched on to avoid overplotting.

* `checkPairwise()` gains an argument `excludeInbred`, which is TRUE by default. This is sensible since the plot shows estimated kappa coefficients, which are well-behaved only for pairs of noninbred individuals.

## Other changes

* **forrel** now requires R version 4.1 and recent versions of **pedtools** and **ribd**. This allowed many simplifications in code and examples.

* Added **scales** as a suggested package.


# forrel 1.3.0

## New features

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
