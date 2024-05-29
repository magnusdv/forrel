#' Pairwise relatedness estimation
#'
#' Estimate the IBD coefficients \eqn{\kappa = (\kappa_0, \kappa_1,
#' \kappa_2)}{(k0, k1, k2)} or the condensed identity coefficients \eqn{\Delta =
#' (\Delta_1, ..., \Delta_9)}{(d1, ..., d9)} between a pair (or several pairs)
#' of pedigree members, using maximum likelihood methods. Estimates of
#' \eqn{\kappa} may be visualised with [showInTriangle()].
#'
#' It should be noted that this procedure estimates the *realised* identity
#' coefficients of each pair, i.e., the actual fractions of the autosomes in
#' each IBD state. These may deviate substantially from the theoretical pedigree
#' coefficients.
#'
#' Maximum likelihood estimation of relatedness coefficients originates with
#' Thompson (1975). Optimisation of \eqn{\kappa} is done in the \eqn{(\kappa_0,
#' \kappa_2)}{(k0, k2)}-plane and restricted to the triangle defined by
#' \deqn{\kappa_0 \ge 0, \kappa_2 \ge 0, \kappa_0 + \kappa_2 \le 1}{k0 >= 0, k2
#' >= 0, k0 + k2 <= 1.}
#'
#' Optimisation of \eqn{\Delta} is done in unit simplex of \eqn{R^8}, using the
#' first 8 coefficients.
#'
#' The implementation optimises the log-likelihood using a projected gradient
#' descent algorithm, combined with a version of Armijo line search.
#'
#' When `param = "kappa"`, the output may be fed directly to [showInTriangle()]
#' for visualisation.
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either a vector with ID labels, or a data frame/matrix with two
#'   columns, each row indicating a pair of individuals. The entries are coerced
#'   to characters, and must match uniquely against the ID labels of `x`. By
#'   default, all pairs of genotyped members of `x` are included.
#' @param markers A vector with names or indices of markers attached to x,
#'   indicating which markers to include. By default, all markers are
#'   used.
#' @param param Either "kappa" (default) or "delta"; indicating which set of
#'   coefficients should be estimated.
#' @param start A probability vector (i.e., with nonnegative entries and sum 1)
#'   of length 3 (if `param = "kappa"`) or 9 (if `param = "delta"`), indicating
#'   the initial value of for the optimisation. By default, `start` is set to
#'   `(1/3, 1/3, 1/3)` if `param = "kappa"` and `(1/9, ..., 1/9)` if `param =
#'   "delta"`.
#' @param tol,beta,sigma Control parameters for the optimisation routine; can
#'   usually be left untouched.
#' @param contourPlot A logical. If TRUE, contours of the log-likelihood
#'   function are plotted overlaying the IBD triangle.
#' @param levels (Only relevant if `contourPlot = TRUE`.) A numeric vector of
#'   levels at which to draw contour lines. If NULL (default), the levels are
#'   chosen automatically.
#' @param verbose A logical.
#'
#' @return An object of class `ibdEst`, which is basically a data frame with
#'   either 6 columns (if `param = "kappa"`) or 12 columns (if `param =
#'   "delta"`). The first three columns are `id1` (label of first individual),
#'   `id2` (label of second individual) and `N` (the number of markers with no
#'   missing alleles). The remaining columns contain the coefficient estimates.
#'
#' @author Magnus Dehli Vigeland
#'
#' @seealso [ibdBootstrap()]
#'
#' @references
#'
#' * E. A. Thompson (1975). _The estimation of pairwise relationships._ Annals
#' of Human Genetics 39.
#'
#' * E. A. Thompson (2000). _Statistical Inference from Genetic Data on
#' Pedigrees._ NSF-CBMS Regional Conference Series in Probability and
#' Statistics. Volume 6.
#'
#' @examples
#'
#' ### Example 1: Siblings
#'
#' # Create pedigree and simulate 100 markers
#' x = nuclearPed(2) |> markerSim(N = 100, alleles = 1:4, seed = 123)
#' x
#'
#' # Estimate kappa (expectation: (0.25, 0.5, 0.25)
#' k = ibdEstimate(x, ids = 3:4)
#' k
#'
#' # Visualise estimate
#' showInTriangle(k, labels = TRUE)
#'
#' # Contour plot of the log-likelihood function
#' ibdEstimate(x, ids = 3:4, contourPlot = TRUE)
#'
#'
#' ### Example 2: Full sib mating
#' y = fullSibMating(1) |>
#'   markerSim(ids = 5:6, N = 1000, alleles = 1:10, seed = 123)
#'
#' # Estimate
#' ibdEstimate(y, param = "delta")
#'
#' # Exact coefficient by `ribd`:
#' ribd::condensedIdentity(y, 5:6, simplify = FALSE)
#'
#' @export
ibdEstimate = function(x, ids = typedMembers(x), param = c("kappa", "delta"),
                       markers = NULL, start = NULL, tol = sqrt(.Machine$double.eps),
                       beta = 0.5, sigma = 0.5, contourPlot = FALSE, levels = NULL,
                       verbose = TRUE) {
  st = Sys.time()
  param = match.arg(param)

  if(!is.null(markers))
    x = selectMarkers(x, markers)

  # Argument `ids` should be either a vector or a matrix-like with 2 columns
  if(is.vector(ids))
    ids = .comb2(ids, vec = TRUE)
  else if(is.data.frame(ids))
    ids = as.matrix(ids)
  if(!is.matrix(ids) || nrow(ids) == 0 || ncol(ids) != 2)
      stop2("`ids` must be either a vector of length at least 2, or a matrix-like with 2 columns")
  mode(ids) = "character"

  # Convert to list of pairs
  pairs = lapply(seq_len(nrow(ids)), function(i) ids[i, ])

  allids = unique.default(unlist(pairs))
  if(!all(allids %in% typedMembers(x)))
    stop2("Untyped pedigree member: ", setdiff(allids, typedMembers(x)))

  # Alleles and frequencies
  alleleData = .prepAlleleData(x, ids = allids)

  # Start point
  if(is.null(start))
    start = switch(param, kappa = rep(1/3, 3), delta = rep(1/9, 9))

  if(verbose) {
    message(sprintf("Estimating '%s' coefficients", param))
    message(sprintf("Initial search value: (%s) ", rst(start, 3)))
    message(sprintf("Pairs of individuals: %d ", length(pairs)))
  }

  # Estimate each pair
  resList = lapply(pairs, function(pair) {
    if(verbose)
      message(sprintf("  %s: ", paste(pair, collapse = " vs. ")), appendLF = FALSE)

    est = .PGD(alleleData[pair], param = param, start = start, tol = tol, beta = beta, sigma = sigma)

    if(verbose)
      message(sprintf("estimate = (%s), iterations = %d", rst(est$estimate, 3), est$iterations))

    est
  })

  coefs = do.call(rbind, lapply(resList, function(r) r$estimate))
  ids = do.call(rbind, lapply(resList, function(r) r$ids))
  N = unlist(lapply(resList, function(r) r$nMarkers))

  res = structure(data.frame(ids, N, coefs),
                  names = c("id1", "id2", "N", if(param == "kappa") paste0("k", 0:2) else paste0("d", 1:9)),
                  class = c("ibdEst", "data.frame"))

  if(contourPlot) {
    if(param == "delta")
      stop2("Contour plot is available only for 'kappa' estimation")
    if(verbose)
      message("Preparing contour plot")
    contoursKappaML(x, allids, levels = levels)
  }

  if(verbose)
    message("Total time: ", ftime(st))

  res
}

#' @export
print.ibdEst = function(x, digits = 5, ...) {
  y = as.data.frame(x)
  y[4:ncol(y)] = round(y[4:ncol(y)], digits = digits)
  print(y, ...)
}

# Note: This makes `as.numeric()` work!
#' @export
as.double.ibdEst = function(x, ...) {
 unlist(as.list(x)[-(1:3)])
}

#' @export
`[.ibdEst` = function(x, i, j) {
  y = as.data.frame(x)
  y[i, j, drop = FALSE]
}


.PGD = function(dat = NULL, param, start, tol = sqrt(.Machine$double.eps),
                beta = 0.5, sigma = 0.5, maxit = 500, x = NULL, ids = NULL, verbose = FALSE) {

  # Names (for output)
  pair = names(dat) %||% c("_1", "_2")

  # Simplify running this function on its own, with `x` and `ids` (for debugging purposes)
  if(is.null(dat))
    dat = .prepAlleleData(x, ids = ids)

  # Remove markers with missing data
  dat = .removeMissing(dat)

  # Coordinate-wise likelihoods: P(G_j | IBD = i)
  wei = .likelihoodWeights(dat, param = param)

  # Log-likelihood function: Input full-dimensional kappa or delta
  loglik = function(p, grad = FALSE) {
    liks = as.numeric(p %*% wei)
    if(!grad)
      return(sum(log(liks)))

    list(loglik = sum(log(liks)), grad = as.numeric(wei %*% (1/liks)))
  }

  # Optimise: Projected gradient descent
  k = 0
  xk = start
  LL = loglik(xk, grad = TRUE)
  if(LL$loglik == -Inf)
    stop2("Initial value is impossible")

  ak = 1

  while(TRUE) {
    ll = LL$loglik
    gr = LL$grad

    if(stationary(xk, gr, tol = tol))
      break

    k = k + 1

    ARMIJO = function(y) {
      LHS = loglik(y, grad = FALSE)
      RHS = ll + sigma * as.numeric(gr %*% (y - xk))
      if(verbose)
        message("Armijo: y = ", rst(y,5), "; LHS = ", rst(LHS,5), "; Diff = ", LHS - RHS)
      if(abs(LHS - RHS) < tol)
        return(NA)
      LHS > RHS
    }

    # Armijo-rule line search (improved)
    y = simplexProject(xk + ak * gr)

    if(is.na(arm <- ARMIJO(y)))
      break

    if(arm) {
      if(verbose) message("Armijo OK - trying to increase step size")
      while(TRUE) {
        ak.try = ak/beta
        if(verbose) message("Increasing ak to ", ak.try)
        y.try = simplexProject(xk + ak.try * gr)
        if(all(abs(y.try - y) < tol) || !isTRUE(ARMIJO(y.try)))
          break
        ak = ak.try
        y = y.try
      }
    }
    else {
      if(verbose) message("Armijo FAIL - decreasing step size")
      while(TRUE) {
        if(ak < 1e-50) {stop2("Precision problems encountered")}
        ak = ak * beta
        if(verbose) message("Decreasing ak to ", ak)
        y = simplexProject(xk + ak * gr)
        if(!isFALSE(ARMIJO(y)))
          break
      }
    }

    newLL = loglik(y, grad = TRUE)

    if(verbose)
      message(sprintf("*** Iteration %d: Point = %s; loglik = %s", k, rst(y, 5), rst(newLL$loglik, 5)))

    # Safety stop
    if(k >= maxit)
      break

    xk = y
    LL = newLL
  }

  if(verbose)
    message("Iterations: ", k)

  #res = data.frame(id1 = pair[1], id2 = pair[2], N = sum(keep), rbind(xk), row.names = NULL)
  #names(res)[-(1:3)] = switch(param, kappa = paste0("k", 0:2), delta = paste0("d", 1:9))

  list(estimate = xk, loglik = LL$loglik, iterations = k, ids = pair, nMarkers = length(dat[[1]][[1]]))
}




# Plot contour lines for the ML function when estimating kappa.
# This function is (optionally) called from within `ibdEstimate()`.
#' @importFrom stats quantile
#' @importFrom graphics contour
contoursKappaML = function(x, ids, peak = NA, levels = NULL) {
  ids = as.character(ids)
  if(length(ids) != 2)
    stop2("Contour plots require `ids` to be a single pair of individuals")

  if(isTRUE(is.na(peak))) {
    peak = ibdEstimate(x, ids, param = "kappa", verbose = FALSE)
    peak = c(peak$k0, peak$k2)
  }

  # Log-likelihood function: Input full-dimensional kappa
  loglik = ibdLoglikFUN(x, ids, input = "kappa")

  n = 51
  k0 = seq(0, 1, length.out = n)
  k2 = seq(0, 1, length.out = n)

  loglikMat = matrix(NA_real_, ncol = n, nrow = n)
  for(i in 1:n) for(j in seq_len(n-i)) {
    kap = c(k0[i], 1 - k0[i] - k2[j], k2[j])
    loglikMat[i,j] = loglik(kap)
  }

  if(is.null(levels)) {
    ll = as.numeric(loglikMat)
    ll = ll[ll > -Inf & !is.na(ll)]

    # nice set of quantiles
    levs = quantile(ll, c(0.10, 0.23, 0.36, 0.49, 0.62, 0.75, 0.88, 0.945, 0.985))

    # round maximally without distorting diffs too much
    d = 6
    while(TRUE) {
      rlevs = unique.default(round(levs, d))
      if(length(rlevs) < length(levs))
        break
      diffchange = (diff(rlevs) - diff(levs))/diff(levs)
      if(any(abs(diffchange) > 0.1))
        break
      d = d - 1
    }

    levels = unique.default(round(levs, d+1))
  }

  ribd::ibdTriangle()
  contour(k0, k2, z = loglikMat, add = TRUE, levels = levels)
  if(!is.null(peak))
    ribd::showInTriangle(peak, new = FALSE)
}

add = function(v, col = 2, pch = 16) points(v[1], v[3], col = col, pch = pch)


# Efficient version directly from allele data, e.g. output from
# `profileSimParametric(...,  returnValue = "alleles")`
# NB: Used in `ibdBootstrap()`
.ibdEstimFromAlleles = function(als, freqList, param, start = NULL) {
  # als: list of 4 vectors with true (not internal) alleles
  # freqList: List like `NorwegianFrequencies`

  # Prepare for estimation
  dat = .prepAlleleData2(als, freqList)

  .PGD(dat, param = param, start = start)$estimate
}
