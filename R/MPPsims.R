#' Missing person power simulations
#'
#' Estimate the exclusion/inclusion power for various selections of available
#' individuals.
#' @param reference A connected `ped` object
#' @param missing The label of the missing person. By default "MP".
#' @param selections A list of pedigree member subsets. In the special case that
#'   all subsets consist of a single individual, `selections` can be given as a
#'   simple vector.
#' @param baseline A logical. If TRUE (default) an *empty* selection, named
#'   "Baseline", is added as the first element of `selection`.
#' @param nProfiles The number of profile simulations for each selection.
#' @param lrSims,thresholdIP Parameters passed onto [missingPersonIP()]
#' @param verbose A logical
#'
#' @return A list with the same length as `selections`. Each entry has elements
#'   `ep` and `ip`, each of which is a list of length `nProfiles`.
#' @export
#'
#' @examples
#' x = nuclearPed(fa = "Gf", mo = "Gm", children = c("Uncle", "Mother"), sex = 1:2)
#' x = addChildren(x, fa = "Father", mo = "Mother", nch = 3, sex = c(1,2,1),
#'                 id = c("S1", "S2", "MP"))
#' x = addSon(x, "Father", id = "HS")
#'
#' # Brother S1 is already genotyped with a marker with 5 alleles
#' m = marker(x, S1 = 1:2, alleles = 1:5)
#' x = setMarkers(x, m)
#'
#' # Alternatives for additional genotyping
#' sel = list("Father", "S2", "HS", "Gm", "Uncle", c("Gm", "Uncle"))
#'
#' plot(x, marker = 1, shaded = sel)
#'
#' # Simulate
#' simData = MPPsims(x, selections = sel, nProfiles = 10, lrSims = 10)
#'
#' # Power plot
#' epSims = lapply(simData, '[[', 'ep')
#' ipSims = lapply(simData, '[[', 'ip')
#' powerPlot(epSims, ipSims, type = 3)
#'
MPPsims = function(reference, missing = "MP", selections, baseline = TRUE,
                   nProfiles = 1, lrSims = 1, thresholdIP = NULL, verbose = TRUE) {

  if(!is.ped(reference))
    stop2("Expecting a connected pedigree as H1")

  if(!hasMarkers(reference))
    stop2("No markers attached to `reference` pedigree")

  if(!is.list(selections))
    selections = as.list(selections)

  powSims = lapply(selections, function(ids) {
    if(verbose) message("Selection: ", toString(ids))

    # Simulate profile for `ids`
    sims = profileSim(reference, ids = ids, N = nProfiles)

    # Compute updated EP and IP for each profile
    ep = lapply(sims, function(y)
      missingPersonEP(y, missing = "MP", verbose = FALSE))

    ip = lapply(sims, function(y)
      missingPersonIP(y, missing = "MP", verbose = FALSE, nsim = lrSims, threshold = thresholdIP))

    list(ep = ep, ip = ip)
  })

  names(powSims) = sapply(selections, toString)

  if(baseline) {
    ep0 = missingPersonEP(reference, missing = missing, verbose = FALSE)
    ip0 = missingPersonIP(reference, missing = missing, nsim = lrSims, threshold = thresholdIP, verbose = FALSE)
    pow0 = list(ep = ep0, ip = ip0)
    powSims = c(list(Baseline = pow0), powSims)
  }

  powSims
}
