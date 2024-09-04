#' LR calculations for paternity and sibship
#'
#' A thin wrapper around [kinshipLR()] for the common scenario of testing a pair
#' of individuals for paternity and/or sibship, against being unrelated.
#'
#' @param x A ped object or a list of such.
#' @param ids A vector of two typed members of `x`. If not given, the typed
#'   members of `x` are selected by default, but note that this gives an error
#'   if the number of such individuals is not 2.
#' @param test The hypotheses to be tested (against 'unrelatedness'). Allowed
#'   values are "pat" (=paternity), "sib" (=full siblings), and "half" (=half
#'   siblings). By default, all three are included.
#'
#' @return A (slightly simplified) `LRresult` object, as described in
#'   [kinshipLR()].
#'
#' @examples
#' # Simulate 100 markers for half siblings
#' x = halfSibPed() |> markerSim(N = 100, ids = 4:5, alleles = 1:3, seed = 1)
#'
#' # Test paternity, full sib, half sib
#' quickLR(x)
#'
#' @export
quickLR = function(x, ids = typedMembers(x), test = c("pat", "sib", "half")) {
  if(length(ids) != 2)
    stop2("Please indicate exactly two individuals. Received: ", ids)

  if(!all(ids %in% typedMembers(x)))
    stop2("Individual is not typed: ", setdiff(ids, typedMembers(x)))

  if(!all(test %in% c("pat", "sib", "half")))
    stop2("Illegal value of `test` (allowed are 'pat', 'sib', 'half'): ", setdiff(test, c("pat", "sib", "half")))

  sex = getSex(x, ids)

  peds = list(
    pat  = nuclearPed(father = if(sex[1] == 1) ids[1] else "_fa_",
                      mother = if(sex[1] == 2) ids[1] else "_mo_",
                      children = ids[2], sex = sex[2]),
    sib  = nuclearPed(father = "_fa_", mother = "_mo_", children = ids, sex = sex),
    half = halfSibPed(sex1 = sex[1], sex2 = sex[2]) |>
      relabel(old = 1:5, new = c("_mo1_", "_fa_", "_mo2_", ids))
  )[test]

  peds$unrel = transferMarkers(from = x, to = singletons(ids, sex = sex))
  n = length(peds)

  res = kinshipLR(peds, ref = n, source = n)
  res$LRtotal = res$LRtotal[-n]
  res$LRperMarker = res$LRperMarker[, -n, drop = FALSE]
  res
}
