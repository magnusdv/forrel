#' LR calculations for paternity and sibship
#'
#' A thin wrapper around [kinshipLR()] for the common scenario of testing a pair
#' of individuals for paternity or sibship against being unrelated.
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector of two typed members of `x`. If not given, the typed
#'   members of `x` are selected by default, but note that this gives an error
#'   if the number of such individuals is not 2.
#' @param test The hypotheses to be tested against unrelatedness. Allowed values
#'   are "pat" (= paternity), "sib" (= full siblings), "half" (= half siblings)
#'   and "cous" (= first cousins). By default, all are included.
#'
#' @return A (slightly simplified) `LRresult` object, as described in
#'   [kinshipLR()].
#'
#' @examples
#' # Simulate 100 markers for half siblings
#' x = halfSibPed() |> markerSim(N = 100, ids = 4:5, alleles = 1:3, seed = 1)
#'
#' # Test for paternity, full sibs, half sibs, 1st cousins
#' quickLR(x)
#'
#' @export
quickLR = function(x, ids = typedMembers(x), test = c("pat", "sib", "half", "cous")) {
  if(length(ids) != 2)
    stop2("Please indicate exactly two individuals. Received: ", ids)

  if(!all(ids %in% typedMembers(x)))
    stop2("Individual is not typed: ", setdiff(ids, typedMembers(x)))

  rels = c("pat", "sib", "half", "cous")
  if(!all(test %in% rels))
    stop2("Illegal value of `test`: ", setdiff(test, rels),
          "\nAllowed inputs are 'pat', 'sib', 'half', 'cous'")

  sex = getSex(x, ids)

  peds = list(
    pat  = nuclearPed(father = if(sex[1] == 1) ids[1] else "_fa_",
                      mother = if(sex[1] == 2) ids[1] else "_mo_",
                      children = ids[2], sex = sex[2]),
    sib  = nuclearPed(father = "_fa_", mother = "_mo_", children = ids, sex = sex),
    half = halfSibPed(sex1 = sex[1], sex2 = sex[2]) |>
      relabel(old = 1:5, new = c("_mo1_", "_fa_", "_mo2_", ids)),
    cous = cousinPed(1) |> setSex(7:8, sex = sex) |>
      relabel(old = 1:8, new = c(paste0("_nn_", 1:6), ids))
  )[test]

  peds$unrel = transferMarkers(from = x, to = singletons(ids, sex = sex))
  n = length(peds)

  res = kinshipLR(peds, ref = n, source = n)
  res$LRtotal = res$LRtotal[-n]
  res$LRperMarker = res$LRperMarker[, -n, drop = FALSE]
  res
}
