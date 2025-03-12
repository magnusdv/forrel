#' Identify and rank the most likely profiles of a pedigree member
#'
#' Markers are assumed independent. For each marker, the possible genotypes for
#' `id` are evaluated and ranked according to how probable they are. It may be
#' wise to try first with `maxPerMarker = 1` to limit computation time,
#' particularly if mutations are modelled.
#'
#' @param A `ped` object with attached markers.
#' @param id Name of the individual to be predicted.
#' @param markers Names or indices of the markers to be included. Default: all.
#' @param maxPerMarker The number of top candidates to be considered per marker.
#'   Default: all.
#' @param verbose A logical, by default TRUE.
#'
#' @return A list with three components. The first, a data frame with profiles
#'   ranked according to likelihood. Under a flat prior, the posterior equals
#'   the likelihood. Then follows the most likely profile. Finally, the second
#'   to most likely genotypes are given.
#'
#' @examples
#'
#' x = nuclearPed(2, father = "FA") |>
#'   addMarker(`3` = "1/1", `4` = "1/2", alleles = 1:2, name = "m1") |>
#'   addMarker(`3` = "1/1", `4` = "2/2", alleles = 1:2, name = "m2") |>
#'   addMarker(`3` = "2/2", `4` = "2/2", alleles = 1:2, name = "m3")
#'
#' rankProfiles(x, "FA")
#' rankProfiles(x, "FA", maxPerMarker = 2)
#' rankProfiles(x, "FA", maxPerMarker = 1)
#'
#' # Same example with mutations allowed
#' y = setMutmod(x, model = "equal", rate = 0.002)
#' rankProfiles(y, "FA")
#'
#' @export
rankProfiles = function(x, id, markers = NULL, maxPerMarker = Inf, verbose = FALSE) {

  if(!hasMarkers(x))
    stop2("No markers attached to pedigree")

  if(is.null(markers))
    markers = seq_len(nMarkers(x))

  # Extract wanted markers
  x = selectMarkers(x, markers)
  nmark = nMarkers(x)

  # Marker names
  mnames = name(x)
  mnames[is.na(mnames)] = as.character((1:nmark)[is.na(mnames)])

  # Likelihood for each marker
  omdList = lapply(1:nmark, function(i) {
    omd = oneMarkerDistribution(x, ids = id, marker = i, verbose = FALSE)
    a = omd[omd > 0, drop = FALSE]
    a = sort.default(a, decreasing = T)

    if(maxPerMarker < length(a))
      a = a[1:maxPerMarker]

    # Fix if only one element
    if(!is.array(a)) a = as.array(a)

    # Trick to ensure proper naming in `expand.grid` below
    names(dimnames(a)) = mnames[i]

    if(verbose) { print(a); cat("\n") }
    a
  })

  # Find most likely profile
  maxGenos = unlist(lapply(omdList, `[`, 1))
  names(maxGenos) = paste(mnames, names(maxGenos), sep=":")

  # Find second most likely marginally
  seconds = lengths(omdList) > 1
  if (any(seconds)){
    marginal2 = unlist(lapply(omdList[seconds], `[`, 2))
    names(marginal2) = paste(mnames[seconds], names(marginal2), sep=":")
  } else
    marginal2 = NA

  # Array with total posterior probs
  pp = Reduce(`%o%`, omdList)

  # Transform dimnames into data frame with possible profiles
  profiles = expand.grid(dimnames(pp))

  # Add probabilities
  prob = as.numeric(pp)
  profiles$likelihood = prob

  # Sort
  profiles = profiles[order(prob, decreasing = T), , drop = F]
  row.names(profiles) = NULL

  list("rankedProfiles" = profiles, "maxProfile" = maxGenos, "no2" = marginal2)
}
