# forrel 1.0.1

## New features

* `readFam()` now has a parameter `Xchrom` which can be used to indicate that
the markers included in the file are on the X chromosome

* `MPPsims()` is more flexible, and allows subsetting of its output.

* `powerPlot()` is more flexible and allows finer control of the plot contents

## Bug fixes

* Fixed several glitches in `readFam()`. It is more robust now, and fails 
gracefully in certain situations which cannot currently be handled (e.g. if 
the file contains twins).


# forrel 1.0.0

* Initial CRAN release
