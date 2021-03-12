#' Pairwise relatedness estimation
#'
#' Estimate the IBD coefficients \eqn{\kappa = (\kappa_0, \kappa_1,
#' \kappa_2)}{(k0, k1, k2)} or the condensed identity coefficients \eqn{\Delta =
#' (\Delta_1, ..., \Delta_9)}{(d1, ..., d9)} between a pair (or several pairs)
#' of pedigree members, using maximum likelihood methods.
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
#' >= 0, k0 + k2 <= 1}. Optimisation of \eqn{\Delta} is done in unit simplex of
#' R^8, using the first 8 coefficients.
#'
#' The implementation optimises the log-likelihood using a projected gradient
#' descent algorithm, combined with a version of Armijo line search.
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either a vector with ID labels, or a data frame/matrix with two
#'   columns, where each row contains the ID labels of two individuals. The
#'   entries are coerced to characters, and must match uniquely against the ID
#'   labels of `x`. By default, all pairs of members of `x` are included.
#' @param markers A vector with names or indices of markers attached to x,
#'   indicating which markers to include. If NULL (default), all markers are
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
#' @return An object of class `ibdEst`, which is basically a dataframe with
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
#' x = nuclearPed(2)
#'
#' # Simulate 100 markers
#' x = markerSim(x, N = 100, alleles = 1:4, seed = 123, verbose = FALSE)
#'
#' # Estimate kappa (expectation: (0.25, 0.5, 0.25)
#' ibdEstimate(x, ids = 3:4)
#'
#' # Plot contours of the log-likelihood function
#' ibdEstimate(x, ids = 3:4, contourPlot = TRUE)
#'
#' ### Example 2: Full sib mating
#' y = fullSibMating(1)
#'
#' # Simulate 200 SNP markers
#' y = markerSim(y, N = 1000, alleles = 1:10, seed = 123, verbose = FALSE)
#'
#' # Estimate
#' ibdEstimate(y, ids = 5:6, param = "delta")
#'
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
  alleleData = .getAlleleData2(x, ids = allids)

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
    message("Total time: ", format(Sys.time() - st, digits = 3))

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

  # Simplify running this function on its own, with `x` and `ids` (for debugging purposes)
  if(is.null(dat))
    dat = .getAlleleData2(x, ids = ids)

  pair = names(dat) %||% c("_1", "_2")

  # Remove missing
  keep = !is.na(dat[[1]]$f1) & !is.na(dat[[2]]$f2)
  if(!all(keep))
    dat = list(lapply(dat[[1]], `[`, keep), lapply(dat[[2]], `[`, keep))

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
      LHS = loglik(y)
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

  list(estimate = xk, loglik = LL$loglik, iterations = k, ids = pair, nMarkers = sum(keep))
}

# Prepare data for fast computation of IBD likelihood
# Input: ids = a vector of ID labels
# Output: List of lists. For each indiv, the output contains a list of 4 vectors:
#   a1: first allele
#   a2: second allele
#   f1: frequency of first allele
#   f2: frequency of second allele
.getAlleleData2 = function(x, ids) {
  nMark = nMarkers(x)
  nSeq = seq_len(nMark)
  ids = as.character(ids)

  if(is.pedList(x)) {
    pednr = getComponent(x, ids, checkUnique = TRUE)
    if(all(pednr == pednr[1]))
      x = x[[pednr[1]]]
  }

  if(is.pedList(x)) {
    alsMat = do.call(rbind, lapply(x, function(cmp) matrix(unlist(cmp$MARKERS), ncol = 2*nMark)))
    freqlist = lapply(nSeq, function(i) unname(afreq(x, i))) # checks consistency!
  }
  else {
    alsMat = matrix(unlist(x$MARKERS), ncol = 2*nMark)
    freqlist = lapply(x$MARKERS, function(m) attr(m, "afreq")) # faster than the generic afreq()
  }

  rownames(alsMat) = unlist(labels(x))
  alsMat = alsMat[ids, , drop = FALSE]

  maxAlNum = max(lengths(freqlist))
  freqMat = matrix(unlist(lapply(freqlist, `length<-`, maxAlNum),
                          recursive = FALSE, use.names = FALSE),
                   ncol = nMark)

  even = 2 * nSeq
  odd = even - 1

  res = lapply(ids, function(id) {
    a1 = alsMat[id, odd]
    a2 = alsMat[id, even]
    a1nonz = a1 > 0
    a2nonz = a2 > 0
    f1 = f2 = rep(NA_real_, nMark)
    f1[a1nonz] = freqMat[cbind(a1, nSeq)]
    f2[a2nonz] = freqMat[cbind(a2, nSeq)]
    list(a1 = a1, a2 = a2, f1 = f1, f2 = f2)
  })

  names(res) = ids
  res
}


# Input: alleleData = ouput from .getAlleleData() for a pair of indivs
.likelihoodWeights = function(alleleData, param = "kappa") {
  .a = alleleData[[1]]$a1
  .b = alleleData[[1]]$a2
  .c = alleleData[[2]]$a1
  .d = alleleData[[2]]$a2
  pa = alleleData[[1]]$f1
  pb = alleleData[[1]]$f2
  pc = alleleData[[2]]$f1
  pd = alleleData[[2]]$f2

  homoz1 = .a == .b
  homoz2 = .c == .d
  mac = .a == .c
  mbc = .b == .c
  mad = .a == .d
  mbd = .b == .d
  g1.fr = 2^(!homoz1) * pa * pb
  g2.fr = 2^(!homoz2) * pc * pd

  if(param == "kappa") {
    # Matrix with 3 rows: conditional probs given UN, PO, MZ
    w = matrix(0, nrow = 3, ncol = length(.a))
    w[1, ] = g1.fr * g2.fr                  # Prob(g1, g2 | UN)
    w[2, ] = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd)) # Prob(g1, g2 | PO)
    w[3, ] = g1.fr * ((mac & mbd) | (mad & mbc)) # Prob(g1, g2 | MZ)
  }
  else if(param == "delta") {
    # Matrix with 9 rows: Row i are the conditional probs given state Ji
    w = matrix(0, nrow = 9, ncol = length(.a))
    w[1, ] = pa * (homoz1 & homoz2 & mac)   # Prob(g1, g2 | J1)
    w[2, ] = pa * pc * (homoz1 & homoz2)    # Prob(g1, g2 | J2)
    w[3, ] = pa * homoz1 * (mac *pd + mad * pc) * (.5)^homoz2 # Prob(g1, g2 | J3)
    w[4, ] = pa * g2.fr * homoz1            # Prob(g1, g2 | J4)
    w[5, ] = pc * homoz2 * (mac *pb + mbc * pa) * (.5)^homoz1 # Prob(g1, g2 | J5)
    w[6, ] = g1.fr * pc * homoz2            # Prob(g1, g2 | J6)
    w[7, ] = g1.fr * ((mac & mbd) | (mad & mbc)) # Prob(g1, g2 | MZ)
    w[8, ] = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd)) # Prob(g1, g2 | PO)
    w[9, ] = g1.fr * g2.fr                  # Prob(g1, g2 | UN)
  }
  else
    stop2("`param` must be either 'kappa' or 'delta': ", param)

  w
}


stationary = function(p, grad, tol) {
  # print(grad,digits = 12)
  # print(sapply(seq_along(p), function(i) {pi = p; pi[i] = p[i] - 1; grad %*% pi}))
  for(i in seq_along(p)) {
    pi = p
    pi[i] = p[i] - 1
    if(grad %*% pi < -tol)
      return(FALSE)
  }
  TRUE
}

# Project any vector in R^D onto the probability simplex.
# Algorithm found in Wang et al. (2013): Projection onto the probability simplex.
simplexProject = function(y) {
  D = length(y)
  if(D == 3)
    u = .sort3(y)
  else
    u = sort.int(y, decreasing = TRUE, method = "shell")
  v = u + 1/seq_len(D) * (1 - cumsum(u))
  p = match(TRUE, v <= 0, nomatch = D + 1) - 1
  lambda = 1/p * (1 - sum(u[1:p]))
  x = y + lambda
  x[x < 0] = 0
  x
}

.sort3 = function(y) {
  if(y[1] < y[2]) {tmp = y[1]; y[1] = y[2]; y[2] = tmp}
  if(y[2] < y[3]) {tmp = y[2]; y[2] = y[3]; y[3] = tmp}
  if(y[1] < y[2]) {tmp = y[1]; y[1] = y[2]; y[2] = tmp}
  y
}




# Plot contour lines for the ML function when estimating kappa.
# This function is (optionally) called from within `ibdEstimate()`.
#' @importFrom stats quantile
contoursKappaML = function(x, ids, peak = NA, levels = NULL) {
  ids = as.character(ids)
  if(length(ids) != 2)
    stop2("Contour plots require `ids` to be a single pair of individuals")

  dat = .getAlleleData2(x, ids = ids)

  if(isTRUE(is.na(peak))) {
    peak = ibdEstimate(x, ids, param = "kappa", verbose = FALSE)
    peak = c(peak$k0, peak$k2)
  }

  # Remove missing
  keep = !is.na(dat[[1]]$f1) & !is.na(dat[[2]]$f2)
  if(!all(keep))
    dat = list(lapply(dat[[1]], `[`, keep), lapply(dat[[2]], `[`, keep))

  # Coordinate-wise likelihoods: P(G_j | IBD = i)
  wei = .likelihoodWeights(dat, param = "kappa")

  # Log-likelihood function: Input full-dimensional kappa
  loglik = function(k0, k2) sum(log(as.numeric(c(k0, 1 - k0 - k2, k2) %*% wei)))

  n = 51
  k0 = seq(0, 1, length.out = n)
  k2 = seq(0, 1, length.out = n)

  loglikMat = matrix(NA_real_, ncol = n, nrow = n)
  for(i in 1:n) for(j in seq_len(n-i))
    loglikMat[i,j] = loglik(k0[i], k2[j])

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

# round + toSTring
rst = function(v, digits = 3) toString(round(v, digits))
