#' Expected likelihood ratio
#'
#' This function computes the expected LR for a single marker, in a kinship test
#' comparing two hypothesised relationships between a set of individuals. The
#' true relationship may differ from both hypotheses. Some individuals may
#' already be genotyped, while others are available for typing. The
#' implementation uses `oneMarkerDistribution()` to find the joint genotype
#' distribution for the available individuals, conditional on the known data, in
#' each pedigree.
#'
#' @param numeratorPed A `ped` object.
#' @param denominatorPed A `ped` object.
#' @param truePed A `ped` object.
#' @param ids A vector of ID labels corresponding to untyped pedigree members.
#'   (These must be members of all three input pedigrees).
#' @param marker either a marker object compatible with `numeratorPed`, or the
#'   name or index of a marker attached to `numeratorPed`.
#'
#' @return A positive number.
#'
#' @examples
#'
#' #---------
#' # Curious example showing that ELR may decrease
#' # by typing additional reference individuals
#' #---------
#'
#' # Numerator ped
#' numPed = nuclearPed(father = "fa", mother = "mo", child = "ch")
#'
#' # Denominator ped: fa, mo, ch are unrelated. (Hack!)
#' denomPed = halfSibPed() |> relabel(old = 1:3, new = c("mo", "fa", "ch"))
#'
#' # Scenario 1: Only mother is typed; genotype 1/2
#' p = 0.9
#' m1 = marker(numPed, mo = "1/2", afreq = c("1" = p, "2" = 1-p))
#' expectedLR(numPed, denomPed, ids = "ch", marker = m1)
#'
#' 1/(8*p*(1-p)) + 1/2 # exact formula
#'
#' # Scenario 2: Include father, with genotype 1/1
#' m2 = m1
#' genotype(m2, id = "fa") = "1/1"
#' expectedLR(numPed, denomPed, ids = "ch", marker = m2)
#'
#' 1/(8*p*(1-p)) + 1/(4*p^2) # exact formula
#'
#' @importFrom pedprobr oneMarkerDistribution
#' @export
expectedLR = function(numeratorPed, denominatorPed, truePed = numeratorPed, ids, marker) {

  if(!is.ped(numeratorPed))
    stop2("Argument `numeratorPed` must be a connected `ped` object")
  if(!is.ped(denominatorPed))
    stop2("Argument `denominatorPed` must be a connected `ped` object")
  if(!is.ped(truePed))
    stop2("Argument `truePed` must be a connected `ped` object")

  # Wrapper (for simpler code)
  OMD = function(ped) oneMarkerDistribution(ped, partialmarker = 1, ids = ids, verbose = FALSE)

  # Numerator
  if(is.marker(marker))
    numeratorPed = setMarkers(numeratorPed, marker)
  else
    numeratorPed = selectMarkers(numeratorPed, marker)
  num = OMD(numeratorPed)

  denominatorPed = transferMarkers(from = numeratorPed,
                                   to = denominatorPed,
                                   erase = TRUE)
  den = OMD(denominatorPed)

  # True pedigree
  if(identical(truePed, numeratorPed))
    true = num
  else if(identical(truePed, denominatorPed))
    true = den
  else {
    truePed = transferMarkers(from = numeratorPed,
                              to = truePed,
                              erase = TRUE)
    true = OMD(truePed)
  }

  ELR = sum(true * num/den)
  ELR
}

