#' Find the most likely profiles of a pedigree member
#'
#' Identify and rank the most likely DNA profiles of a pedigree member. For each
#' marker, the possible genotypes of the indicated person are ranked by
#' likelihood.
#'
#' Note that this function assumes that all markers are independent.
#'
#' If the marker data includes mutation models, it may be wise to try first with
#' `maxPerMarker = 1` to limit computation time.
#'
#' @param x A `ped` object with attached markers.
#' @param id The name of a single (typically untyped) pedigree member.
#' @param markers Names or indices of the markers to be included. Default: all.
#' @param maxPerMarker A single number, limiting the number of top genotypes
#'   considered for each marker. Default: `Inf` (no restriction).
#' @param verbose A logical, by default FALSE.
#'
#' @return A list with the following components (`N` denotes the number of
#'   markers):
#'
#' * `profiles`: A data frame with `N+1` columns, containing the possible
#'   profiles, ranked by likelihood.
#' * `marginal1`: A numeric of length `N`, giving the marginal
#'   probability of the most likely genotype for each marker.
#' * `marginal2`: A numeric of length `N`, with marginals for the *second* most
#'   likely genotype for each marker, or `NA` if there is no second.
#' * `best`: A character of length `N` containing the most likely profile.
#'   This is the same as `names(marginal1)`, and also as `profiles[1, 1:N]`.
#'
#' @examples
#'
#' x = nuclearPed(nch = 4) |>
#'    markerSim(N = 4, alleles = c("a", "b", "c"), seed = 1, verbose = FALSE)
#' x
#' # Remove data for father
#' y = setAlleles(x, ids = 1, alleles = 0)
#'
#' # Most likely profiles of father
#' rankProfiles(y, id = 1)
#'
#' # Compare with truth
#' getGenotypes(x, ids = 1)
#'
#' # Same example with mutations allowed
#' z = setMutmod(y, model = "equal", rate = 0.01)
#' rankProfiles(z, id = 1)
#'
#' @export
rankProfiles = function(x, id, markers = NULL, maxPerMarker = Inf, verbose = FALSE) {

  if(!hasMarkers(x))
    stop2("No markers attached to pedigree")

  if(length(id) != 1)
    stop2("Please indicate a single individual")

  if(!is.null(markers))
    x = selectMarkers(x, markers)

  nmark = nMarkers(x)
  if(nmark == 0)
    stop2("No markers selected")

  # Marker names
  mnames = name(x)
  mnames[is.na(mnames)] = as.character((1:nmark)[is.na(mnames)])

  # Likelihood for each marker
  omdList = lapply(1:nmark, function(i) {
    omd = oneMarkerDistribution(x, ids = id, marker = i, verbose = FALSE)
    a = omd[omd > 0, drop = FALSE]
    a = sort.default(a, decreasing = T)

    if(maxPerMarker < length(a))
      length(a) = maxPerMarker

    # Fix if only one element
    if(!is.array(a)) a = as.array(a)

    # Trick to ensure proper naming in `expand.grid` below
    names(dimnames(a)) = mnames[i]

    if(verbose) { print(a); cat("\n") }
    a
  })

  # Most likely profile
  marginal1 = unlist(lapply(omdList, `[`, 1))
  best = names(marginal1)

  # Find second most likely marginally (or NAs)
  marginal2 = unlist(lapply(omdList, `[`, 2))

  # Array with total posterior probs
  pp = Reduce(`%o%`, omdList)

  # Transform dimnames into data frame with possible profiles
  profiles = expand.grid(dimnames(pp), stringsAsFactors = FALSE)

  # Add probabilities
  prob = as.numeric(pp)
  profiles$likelihood = prob

  # Sort
  profiles = profiles[order(prob, decreasing = T), , drop = F]
  row.names(profiles) = NULL

  list(profiles = profiles, best = best, marginal1 = marginal1, marginal2 = marginal2)
}
