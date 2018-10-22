#' Simulation of DNA profiles
#'
#' Simulation of DNA profiles, given a list of marker loci. Some pedigree
#' members may already be genotyped; in that case the simulation is conditional
#' on these.
#'
#' @param x A `ped` object
#' @param N The number of profiles to be simulated
#' @param ids A character (or coercible to character) with ID labels indicating
#'   whose genotypes should be simulated
#' @param conditions A list of marker objects, (or an integer vector indicating
#'   markers attached to `x`). If these contain genotypes the simulations will
#'   condition on these. Locus annotations (allele frequencies, mutationmodels
#'   a.s.o.) are extracted from each marker.
#' @param ... Further arguments passed on to [markerSim()]
#'
#' @return A list of length `N`. Each element is a `ped` object with `K`
#'   attached markers, where `K = length(condition)`.
#'
#' @examples
#' # Example with two brothers, one of which is already genotyped with 2 markers.
#' x = nuclearPed(children = c("B1", "B2"))
#'
#' m1 = marker(x, B1 = 1:2, alleles = 1:3, afreq = c(.2, .3, .5))
#' m2 = marker(x, B1 = 4, alleles = 1:4, afreq = c(.1,.2,.3,.4))#, chrom="X")
#'
#' # These contain the profile of B1
#' cond = list(m1, m2)
#'
#' # Simulate 3 profiles of B2 conditional on the above
#' profileSim(x, N = 3, ids = "B2", condition = cond)
#'
#'
#' @export
profileSim = function(x, N = 1, ids = NULL, conditions = NULL, ...){

  # Iterate over the loci, make N simulations of each.
  sims_markerwise = lapply(conditions, function(pm)
      markerSim(x, N = N, ids = ids, partialmarker = pm, verbose=F, ...))

  # Transpose: Extract i'th marker from each sim above.
  # Output: List of length N, each with length(conditions) markers
  sims = lapply(1:N, function(i) {
      mlist = lapply(sims_markerwise, function(y) y$markerdata[[i]])
      setMarkers(x, mlist)
    })

  sims
}