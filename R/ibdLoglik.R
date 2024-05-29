#' Pairwise IBD likelihood
#'
#' Given genotype data from two individuals, computes the log-likelihood of a
#' single set of IBD coefficients, either  `kappa` = \eqn{(\kappa_0, \kappa_1,
#' \kappa_2)} or the Jacquard coefficients `delta` = \eqn{(\Delta_1, ..., \Delta_9)}.
#' The `ibdLoglikFUN` version returns an efficient _function_ for computing such
#' likelihoods, suitable for optimisations such as in [ibdEstimate()].
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector of ID labels.
#' @param kappa A probability vector of length 3.
#' @param delta A probability vector of length 9.
#' @param input Either "kappa", "kappa02" or "delta". See Value.

#' @return `ibdLoglik()` returns a single number; the total log-likelihood over
#'   all markers included.
#'
#'   `ibdLoglikFUN()` returns a function for computing such log-likelihoods. The
#'   function takes a single input vector `p`, whose interpretation depends on
#'   the `input` parameter:
#'
#'  * "kappa": `p` is expected to be a set of kappa coefficients
#'   \eqn{(\kappa_0, \kappa_1, \kappa_2)}.
#'  * "kappa02": `p` should be a vector of length 2 containing the coefficients
#'   \eqn{\kappa_0} and \eqn{\kappa_2}. This is sometimes a convenient shortcut
#'   when working in the IBD triangle.
#'  * "delta": Expects `p` to be a set of condensed Jacquard coefficients
#'   \eqn{(\Delta_1, ..., \Delta_9)}.
#'
#' @examples
#' # Siblings typed with 10 markers
#' x = nuclearPed(2) |> markerSim(N = 10, alleles = 1:4)
#'
#' # Calculate log-likelihood at a single point
#' k = c(0.25, 0.5, 0.25)
#' ibdLoglik(x, ids = 3:4, kappa = k)
#'
#' # Or first get a function, and then apply it
#' llFun = ibdLoglikFUN(x, ids = 3:4, input = "kappa")
#' llFun(k)
#'
#'
#' @export
ibdLoglik = function(x = NULL, ids = NULL, kappa = NULL, delta = NULL) {
  if(!is.null(kappa)) {
    if(length(kappa) !=  3)
      stop2("Argument `kappa` must have length 3")
    p = kappa
    input = "kappa"
  }
  else if((!is.null(delta))) {
    if(length(delta) !=  9)
      stop2("Argument `delta` must have length 9")
    p = delta
    input = "delta"
  }
  else
    stop2("Either `kappa` or `delta` must be given")

  if(!isTRUE(all.equal(sum(p), 1)))
    stop2("Coefficients do not sum to 1: ", p)

  llFun = ibdLoglikFUN(x = x, ids = ids, input = input)
  llFun(p)
}

#' @rdname ibdLoglik
#' @export
ibdLoglikFUN = function(x, ids, input = c("kappa", "kappa02", "delta")) {

  if(length(ids) != 2)
    stop2("`ids` must have length 2")

  input = match.arg(input)
  dat = .prepAlleleData(x, ids = ids) |>
    .removeMissing()

  # Coordinate-wise likelihoods: P(G_j | IBD = i)
  param = if(input == "kappa02") "kappa" else input
  W = .likelihoodWeights(dat, param = param)

  # Define function
  if(input == "kappa02")
    loglik = function(p) sum(log(as.numeric(c(p[1], 1-p[1]-p[2], p[2]) %*% W)))
  else
    loglik = function(p) sum(log(as.numeric(p %*% W)))

  loglik
}


# Efficient versions directly from alleles, e.g. output of profileSimParametric
# Input: Similar to `.ibdEstimFromAlleles()`

.ibdLoglikFromAlleles = function(als, freqList, kappa = NULL, delta = NULL) {
  if(!is.null(kappa)) {
    if(length(kappa) !=  3)
      stop2("Argument `kappa` must have length 3")
    p = kappa
    input = "kappa"
  }
  else if((!is.null(delta))) {
    if(length(delta) !=  9)
      stop2("Argument `delta` must have length 9")
    p = delta
    input = "delta"
  }
  else
    stop2("Either `kappa` or `delta` must be given")

  if(!isTRUE(all.equal(sum(p), 1)))
    stop2("Coefficients do not sum to 1: ", p)

  llFun = .ibdLoglikFromAllelesFUN(als, freqList, input = input)
  llFun(p)
}


.ibdLoglikFromAllelesFUN = function(als, freqList, input = c("kappa", "kappa02", "delta")) {

  input = match.arg(input)
  dat = .prepAlleleData2(als = als, freqList = freqList) |>
    .removeMissing()

  # Coordinate-wise likelihoods: P(G_j | IBD = i)
  param = if(input == "kappa02") "kappa" else input
  W = .likelihoodWeights(dat, param = param)

  # Define function
  if(input == "kappa02")
    loglik = function(p) sum(log(as.numeric(c(p[1], 1-p[1]-p[2], p[2]) %*% W)))
  else
    loglik = function(p) sum(log(as.numeric(p %*% W)))

  loglik
}
