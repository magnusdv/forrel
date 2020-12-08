#' Random person exclusion power
#'
#' This is a special case of [exclusionPower()], computing the power to exclude
#' a random person as a given pedigree member. More specifically, the function
#' computes the probability of observing, in an individual unrelated to the
#' family individual, a genotype incompatible with the typed family members.
#'
#' @param x A `ped` object with attached markers.
#' @param id The ID label of a single pedigree member.
#' @param disableMutations This parameter determines how mutation models are
#'   treated. Possible values are as follows:
#'
#'   * `NA` (the default): Mutations are disabled only for those markers whose
#'   known genotypes are consistent with the pedigree. This is determined by
#'   temporarily removing all mutation models and checking which markers have
#'   nonzero likelihood.
#'
#'   * `TRUE`: Mutations are disabled for all markers. This will result in an
#'   error if any markers are inconsistent.
#'
#'   * `FALSE`: No action is done to disable mutations.
#'
#'   * A vector containing the names or indices of those markers for which
#'   mutations should be disabled.
#' @inheritParams exclusionPower
#'
#' @return The `EPresult` object returned by [exclusionPower()].
#'
#' @examples
#'
#' # Four siblings:
#' x = nuclearPed(4)
#'
#' # First 3 sibs typed with 4 triallelic markers
#' x = markerSim(x, N = 4, ids = 3:5, alleles = 1:3, seed = 577, verbose = FALSE)
#'
#' # Probability that a random man is excluded as the fourth sibling
#' randomPersonEP(x, id = 6)
#'
#' @importFrom pedprobr likelihood
#' @export
randomPersonEP = function(x, id, markers = NULL, disableMutations = NA, verbose = TRUE) {

  if(!is.ped(x))
    stop2("Expecting a connected pedigree")

  claimPed = relabel(x, old = id, new = "RAN")
  truePed = list(x, singleton("RAN", sex = getSex(x, id)))

  ep = exclusionPower(claimPed = claimPed, truePed = truePed, ids = "RAN",
                      markers = markers, source = "claim",
                      disableMutations = disableMutations,
                      plot = FALSE, verbose = verbose)

  # Change the `ids` entry from "RAN" to `id`
  ep$params$ids = id

  ep
}

