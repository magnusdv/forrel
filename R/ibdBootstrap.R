#' Bootstrap estimation of IBD coefficients
#'
#' This function produces (parametric or nonparametric) bootstrap estimates of
#' the IBD coefficients between two individuals. Both kappa and delta
#' coefficients are supported (see [ibdEstimate()]).
#'
#' The parameter `method` controls how bootstrap estimates are obtained in each
#' replication.
#'
#' If `method = "parametric"`, new profiles for two individuals are simulated
#' from the input coefficients, followed by a re-estimation of the coefficients.
#'
#' If `method = "nonparametric"`, the original markers are sampled with
#' replacement, before the coefficients are re-estimated.
#'
#' @param x A `ped` object. If `method = "parametric"`, this is only used to
#'   extract the allele frequencies, and can be skipped if `freqList` is
#'   provided.
#' @param ids A pair of ID labels.
#' @param param Either NULL (default), "kappa" or "delta". (See below.)
#' @param kappa,delta Probability vectors of length 3 (kappa) or 9 (delta).
#'   Exactly one of `param`, `kappa` and `delta` must be non-NULL. If `kappa`
#'   and `delta` are both NULL, the appropriate set of coefficients is computed
#'   as `ibdEstimate(x, ids, param)`.
#' @param method Either "parametric" or "nonparametric" (abbreviations are
#'   allowed). see Details.
#' @param N The number of simulations.
#' @param freqList A list of probability vectors: The allele frequencies for
#'   each marker.
#' @param plot A logical, only relevant for bootstraps of kappa. If TRUE, the
#'   bootstrap estimates are plotted in the IBD triangle.
#' @param seed An integer seed for the random number generator (optional).
#'
#' @return A data frame with `N` rows containing the bootstrap estimates. The
#'   last column (`dist`) gives the euclidean distance to the original
#'   coefficients, viewed as a point in R^3 (kappa) or R^9 (delta).
#'
#' @seealso [ibdEstimate()]
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
#' boot1 = ibdBootstrap(kappa = c(0.25, .5, .25), N = N, freqList = freqList)
#' boot1
#'
#' # Mean deviation
#' mean(boot1$dist)
#'
#' # Same, but with the 9 identity coefficients.
#' delta = c(0, 0, 0, 0, 0, 0, .25, .5, .25)
#' boot2 = ibdBootstrap(delta = delta, N = N, freqList = freqList)
#'
#' # Mean deviation
#' mean(boot2$dist)
#'
#' #### Non-parametric bootstrap.
#' # Requires `x` and `ids` to be provided
#'
#' x = nuclearPed(2)
#' x = markerSim(x, ids = 3:4, N = 50, alleles = 1:10, seed = 123)
#'
#' bootNP = ibdBootstrap(x, ids = 3:4, param = "kappa", method = "non", N = N)
#'
#' # Parametric bootstrap can also be done with this syntax
#' bootP = ibdBootstrap(x, ids = 3:4, param = "kappa", method = "par", N = N)
#'
#' @export
ibdBootstrap = function(x = NULL, ids = NULL, param = NULL, kappa = NULL, delta = NULL,
                        N, method = "parametric", freqList = NULL, plot = TRUE, seed = NULL) {

  if(sum(!is.null(param), !is.null(kappa), !is.null(delta)) != 1)
    stop2("Exactly one of `param`, `kappa` and `delta` must be given")

  if(!is.null(ids) && length(ids) != 2)
    stop2("`ids` must be either NULL or a vector of length 2")

  if(!is.null(param)) {
    if(is.null(x))
      stop2("The argument `param` requires that `x` and `ids` are also provided")

    if(param == "kappa")
      kappa = ibdEstimate(x, ids = ids, param = "kappa", verbose = FALSE)
    else if(param == "delta")
      delta = ibdEstimate(x, ids = ids, param = "delta", verbose = FALSE)
  }
  else {
    param = if(!is.null(kappa)) "kappa" else "delta"
  }

  # Allow inputs directly from `ibdEstimate()`
  if(inherits(kappa, "ibdEst"))
    kappa = as.numeric(kappa)
  if(inherits(delta, "ibdEst"))
    delta = as.numeric(delta)

  coefs = kappa %||% delta
  if(!isTRUE(all.equal(sum(coefs), 1)))
    stop2("Coefficients do not sum to 1: ", coefs)

  coefNms = if(param == "kappa") paste0("k", 0:2) else paste0("d", 1:9)

  method = match.arg(method, c("parametric", "nonparametric"))

  if(!is.null(seed))
    set.seed(seed)

  if(method == "parametric") {
    freqList = freqList %||% getFreqDatabase(x)

    # Insert default allele labels (1,2,3,...) where needed
    noNames = vapply(freqList, function(fr) is.null(names(fr)), logical(1))
    freqList[noNames] = lapply(freqList[noNames], function(fr) {names(fr) = seq_along(fr); fr})

    # Simulate genotypes
    sims = profileSimParametric(kappa = kappa, delta = delta, N = N, freqList = freqList,
                                returnValue = "internal")

    # Bootstrap estimates
    boots = lapply(sims, function(s)
      .ibdEstimFromAlleles(s, freqList, param = param, start = coefs))
  }
  else {
    if(is.null(x))
      stop2("For non-parametric bootstrap, `x` and `ids` must be provided")

    nMark = nMarkers(x)

    boots = replicate(N, simplify = FALSE, {
      markers = sample.int(nMark, replace = TRUE)
      bs = ibdEstimate(x, ids = ids, param = param, markers = markers, verbose = FALSE)
      as.numeric(bs)
    })
  }

  # Wide matrix - makes "dist" calulation below simpler
  resWide = matrix(unlist(boots, use.names = FALSE), ncol = N,  # faster than do.call(cbind, boots)
                   dimnames = list(coefNms, NULL))

  # Transpose and convert to data.frame
  res = as.data.frame(t.default(resWide))

  # Add column with Euclidean distance from given coefs
  res$dist = sqrt(colSums((resWide - coefs)^2))

  # Plot
  if(plot && param == "kappa") {
    showInTriangle(res, lwd = 1, pch = 1, col = 4)
    showInTriangle(kappa, new = FALSE, col = "red", pch = 4, lwd = 4, cex = 2.2)
  }

  res

}

.ibdEstimFromAlleles = function(als, freqList, param, start = NULL) {
  # NB: als is a list of 4 vectors with true (not internal) alleles
  # Therefore: Need the allele names in each freqList vector

  alsMat = do.call(rbind, als)
  freqMat = vapply(seq_along(freqList),
                   function(i) freqList[[i]][alsMat[, i]],
                   FUN.VALUE = numeric(4))

  # Same format as .getAlleleData2()
  dat = list(list(a1 = als$a, a2 = als$b, f1 = freqMat[1, ], f2 = freqMat[2, ]),
             list(a1 = als$c, a2 = als$d, f1 = freqMat[3, ], f2 = freqMat[4, ]))

  .PGD(dat, param = param, start = start)$estimate
}


#' @export
nonParamBoot = function(x, ids, param = "kappa", N, plot = TRUE, seed = NULL) {

  coefs = as.numeric(ibdEstimate(x, ids = ids, param = param, verbose = FALSE))
  coefNms = names(coefs)

  nMark = nMarkers(x)

  boots = replicate(N, simplify = FALSE, {
    markers = sample.int(nMark, replace = TRUE)
    bs = ibdEstimate(x, ids = ids, param = param, markers = markers, verbose = FALSE)
    as.numeric(bs)
  })

  # Wide matrix - makes "dist" calculation below simpler
  resWide = matrix(unlist(boots, use.names = FALSE),
                   ncol = N, dimnames = list(coefNms, NULL))

  # Transpose and convert to data.frame
  res = as.data.frame(t.default(resWide))

  # Add column with Euclidean distance from given coefs
  res$dist = sqrt(colSums((resWide - coefs)^2))

  # Plot
  if(plot && param == "kappa") {
    showInTriangle(res, lwd = 1, pch = 1, col = 4)
    showInTriangle(coefs, new = FALSE, col = "red", pch = 4, lwd = 4, cex = 2.2)
  }

  res
}
