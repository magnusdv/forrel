#' Bootstrap estimation of IBD coefficients
#'
#' These functions produce (parametric) bootstrap estimates of the IBD
#' coefficients between two individuals. Both kappa coefficients (for noninbred
#' individuals) and the 9 condensed identity coefficients are supported.
#'
#' In each replication, profiles for two individuals are simulated under the
#' input coefficients, and their relationship is re-estimated with [ibdEstim()].
#'
#'
#' @param kappa,delta A probability vector of length 3 (kappa) or 9 (delta): The
#'   coefficients under which simulations are performed.
#' @param N The number of simulations.
#' @param freqList A list of probability vectors: The allele frequencies for
#'   each marker.
#' @param plot A logical. If TRUE, the bootstrap kappa estimates are plotted in
#'   the IBD triangle.
#'
#' @return A data frame with `N` rows containing the bootstrap estimates. The
#'   last column (`dist`) gives the euclidean distance to the input, viewed as a
#'   point in R^3 (kappa) or R^9 (delta).
#'
#' @examples
#'
#' # Frequency list of 15 standard STR markers
#' freqList = NorwegianFrequencies[1:15]
#'
#' # Number of bootstrap simulations (increase!)
#' N = 5
#'
#' # Bootstrap estimates for kappa of full siblings
#' boot1 = kappaBootstrap(c(0.25, .5, .25), N = N, freqList = freqList)
#' boot1
#'
#' # Mean deviation
#' mean(boot1$dist)
#'
#' # Same, but with the 9 identity coefficients.
#' delta = c(0, 0, 0, 0, 0, 0, .25, .5, .25)
#' boot2 = deltaBootstrap(delta, N = N, freqList = freqList)
#'
#' # Mean deviation
#' mean(boot2$dist)
#'
#' @export
kappaBootstrap = function(kappa, N, freqList, plot = TRUE) {

  # Insert default allele labels (1,2,3,...) where needed
  noNames = vapply(freqList, function(fr) is.null(names(fr)), logical(1))
  freqList[noNames] = lapply(freqList[noNames], function(fr) {names(fr) = seq_along(fr); fr})

  # Simulate genotypes
  sims = profileSimParametric(kappa = kappa, N = N, freqList = freqList, returnValue = "alleles")

  # Bootstrap estimates
  boots = do.call(rbind, lapply(sims, function(als)
    .ibdEstimFromAlleles(als, freqList, param = "kappa", start = kappa)))

  # Plot
  if(plot) showInTriangle(boots)

  coeffCols = boots[paste0("k", 0:2)]

  # Euclidean distance from given coeffs
  dist = sqrt(colSums((t.default(as.matrix(coeffCols)) - kappa)^2))

  # Build output
  res = cbind(coeffCols, dist = dist)

  # Keep convergence issues, if any
  if(any(boots$convergence != "OK"))
    res$conv = boots$convergence

  res
}

#' @rdname kappaBootstrap
#' @export
deltaBootstrap = function(delta, N, freqList) {

  # Insert default allele labels (1,2,3,...) where needed
  noNames = vapply(freqList, function(fr) is.null(names(fr)), logical(1))
  freqList[noNames] = lapply(freqList[noNames], function(fr) {names(fr) = seq_along(fr); fr})

  # Simulate genotypes
  sims = profileSimParametric(delta = delta, N = N, freqList = freqList, returnValue = "alleles")

  # Bootstrap estimates
  boots = do.call(rbind, lapply(sims, function(als)
    .ibdEstimFromAlleles(als, freqList, param = "delta", start = delta)))

  coeffCols = boots[paste0("d", 1:9)]

  # Euclidean distance from given coeffs
  dist = sqrt(colSums((t.default(as.matrix(coeffCols)) - delta)^2))

  # Build output
  res = cbind(coeffCols, dist = dist)

  # Keep convergence issues, if any
  if(any(boots$convergence != "OK"))
    res$conv = boots$convergence

  res
}

.ibdEstimFromAlleles = function(als, freqList, param, start = NULL) {
  # NB: als is a list of 4 vectors with true (not internal) alleles
  # Therefore: Need the allele names in each freqList vector

  param = match.arg(param, c("kappa", "delta"))

  alsMat = do.call(rbind, als)
  freqMat = vapply(seq_along(freqList),
                   function(i) freqList[[i]][alsMat[, i]],
                   FUN.VALUE = numeric(4))

  # Same format as .getAlleleData2()
  dat = list(list(a1 = als$a, a2 = als$b, f1 = freqMat[1, ], f2 = freqMat[2, ]),
             list(a1 = als$c, a2 = als$d, f1 = freqMat[3, ], f2 = freqMat[4, ]))

  if(param == "kappa")
    .kappaEstim(dat, start = start, reltol = 1e-12)
  else
    .deltaEstim(dat, start = start, reltol = 1e-12)
}

