#' Expected likelihood ratio
#'
#' This function computes the expected LR for a single marker, in a kinship test comparing
#' two hypothesised relationships between a set of individuals. The true relationship may
#' differ from both hypotheses. Some individuals may already be genotyped, while others
#' are available for typing. The implementation uses `oneMarkerDistribution()` to find the
#' joint genotype distribution for the available individuals, conditional on the known
#' data, in each pedigree.
#'
#' @param numeratorPed A `ped` object or a list of such.
#' @param denominatorPed A `ped` object or a list of such.
#' @param truePed A `ped` object or a list of such.
#' @param ids A vector of ID labels corresponding to untyped pedigree members. (These must
#'   be members of all three input pedigrees).
#' @param marker The name or index of a marker attached to `numeratorPed`. Alternatively,
#'   a standalone `marker` object compatible with `numeratorPed`.
#'
#' @return A positive number; the expected LR.
#'
#' @examples
#'
#' #---------
#' # Curious example showing that ELR may decrease
#' # by typing additional reference individuals
#' #---------
#'
#' # Numerator ped
#' num = nuclearPed(father = "fa", mother = "mo", child = "ch")
#'
#' # Denominator ped: fa, mo, ch are unrelated
#' den = singletons(c("fa", "mo", "ch"), sex = c(1, 2, 1))
#'
#' # Scenario 1: Only mother is typed; genotype 1/2
#' p = 0.9
#' m = marker(num, mo = "1/2", afreq = c("1" = p, "2" = 1-p))
#' expectedLR(num, den, ids = "ch", marker = m)
#'
#' 1/(8*p*(1-p)) + 1/2  # exact
#'
#' # Scenario 2: Include father, with genotype 1/1
#' genotype(m, id = "fa") = "1/1"
#' expectedLR(num, den, ids = "ch", marker = m)
#'
#' 1/(8*p*(1-p)) + 1/(4*p^2)  # exact
#'
#' @importFrom pedprobr oneMarkerDistribution
#' @export
expectedLR = function(numeratorPed, denominatorPed, truePed = numeratorPed, ids, marker) {

  # Wrapper (for simpler code)
  OMD = function(ped) oneMarkerDistribution(ped, marker = 1, ids = ids, verbose = FALSE)

  # Numerator
  if(is.marker(marker)) {
    if(!is.ped(numeratorPed))
      stop2("When `marker` is a standalone marker object, `numeratorPed` must be connected")
    numeratorPed = setMarkers(numeratorPed, marker)
  }
  else
    numeratorPed = selectMarkers(numeratorPed, marker)
  num = OMD(numeratorPed)

  # Check sex if X-linked
  if(isXmarker(numeratorPed)) {
    numSex = getSex(numeratorPed, ids)
    denomSex = getSex(denominatorPed, ids)
    trueSex = getSex(truePed, ids)
    if(!identical(numSex, denomSex) || !identical(numSex, trueSex))
      stop2("For an X-chromosomal marker, the sex of `ids` must agree between all hypotheses")
  }

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

  # Expected LR (avoid  0/0)
  keep = true > 0
  ELR = sum(true[keep] * num[keep]/den[keep])

  ELR
}

