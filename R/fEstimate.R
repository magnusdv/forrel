#' Estimation of inbreeding coefficients
#'
#' Estimate the inbreeding coefficient \eqn{f} of pedigree members from their
#' genotypes. The default method uses maximum likelihood, but moment estimators
#' are also available.
#'
#' The ML estimator assumes independent markers and known allele frequencies.
#' For a single marker with allele frequencies \eqn{p_a}, the genotype
#' likelihood \eqn{P(G \mid f)} is:
#'
#' \deqn{P(G = a/a \mid f) = f p_a + (1-f)p_a^2,} \deqn{P(G = a/b \mid f) =
#' 2(1-f)p_a p_b,\quad a \neq b.}
#'
#' For multiple markers, the full likelihood is the product over markers. The
#' function uses `optimise()` to maximise the log-likelihood over \eqn{[0,1]}.
#'
#' The `"simple"` method is a method-of-moments estimator based on observed and
#' expected homozygosity. Let \eqn{L} be the number of markers, \eqn{N_H} the
#' number of observed homozygous genotypes, and \eqn{H^* = \sum_1^L \sum_a
#' p_a^2} the expected number under Hardy--Weinberg equilibrium. The estimator
#' is then:
#'
#' \deqn{\hat{f} = \dfrac{N_H - H^*}{L - H^*}.}
#'
#' The `"ritland"` method implements the estimator of Ritland (1996). For marker
#' \eqn{j}, let \eqn{n_j} be the number of alleles, and define \eqn{T_j} to be
#' \eqn{1/p_a} if the genotype is \eqn{a/a}, and 0 otherwise. The estimator is
#' then:
#'
#' \deqn{\hat{f} = \dfrac{\sum_{j=1}^L (T_j - 1)}{\sum_{j=1}^L (n_j - 1)}
#'       = \dfrac{\sum_{\text{hom}} 1/p_a - L}{A - L},}
#'
#' where \eqn{A = \sum n_j} is the total number of alleles across markers. Note
#' that both this and the "simple" estimator may produce estimates outside the
#' interval \eqn{[0,1]}.
#'
#' The MLE approach is described explicitly in Section 3.2 of Thompson (1986).
#' Ritland (1996) introduces the Ritland estimator, and also mentions both the
#' simple estimator and the MLE. Note that Ritlands formulation of the MLE is
#' slightly different (but equivalent) to Thompson's.
#'
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector with ID labels, by default, all genotyped members of `x`.
#' @param method A character, either "mle" (maximum-likelihood estimate;
#'   default), "simple" (simple method-of-moments estimator), or "ritland" (see
#'   Details).
#' @param markers A vector with names or indices of markers attached to x,
#'   indicating which markers to include. By default, all markers are included.
#'   Missing genotypes are removed before estimation.
#' @param verbose A logical.
#'
#' @return A named numeric vector.
#'
#' @seealso [ibdEstimate()]
#'
#' @references
#'
#' * Thompson, E. A. (1986). Pedigree Analysis in Human Genetics.
#' _Johns Hopkins University Press_.
#' * Ritland, K. (1996). Estimators for pairwise relatedness and individual
#' inbreeding coefficients. _Genetical Research_.
#'
#' @examples
#'
#' # Simulate DNA profile for child of full sibs
#' x = nuclearPed(2, sex = 1:2) |> addSon(3:4) |>
#'   profileSim(ids = 5, markers = NorwegianFrequencies, seed = 123)
#' x
#'
#' # Estimate inbreeding coefficient (expected: 0.0625)
#' fEstimate(x)
#'
#' #-----------------------------
#' # Compare different estimators
#' #-----------------------------
#'
#' # Simulate 200 DNA profiles for a child of half-sibs
#' x = halfSibPed(sex2 = 2) |> addSon(4:5) |>
#'   profileSim(N = 200, ids = 6, markers = NorwegianFrequencies, seed = 123)
#'
#' # Estimates using all three methods
#' fhat = data.frame(
#'   mle = sapply(x, fEstimate, method = "mle"),
#'   sim = sapply(x, fEstimate, method = "simple"),
#'   rit = sapply(x, fEstimate, method = "ritland")
#' )
#'
#' # Mean estimates
#' apply(fhat, 2, mean)
#'
#' # RMSE
#' apply(fhat, 2, function(v) sqrt(mean((v - 0.125)^2)))
#'
#' @export
fEstimate = function(x, ids = typedMembers(x), method = c("mle", "simple", "ritland"),
                     markers = NULL, verbose = FALSE) {
  st = Sys.time()
  method = match.arg(method)

  if(!is.null(markers))
    x = selectMarkers(x, markers)

  ids = as.character(ids)
  if(!all(ids %in% typedMembers(x)))
    stop2("Untyped pedigree member: ", setdiff(ids, typedMembers(x)))

  # Alleles and frequencies
  alleleData = .prepAlleleData(x, ids = ids)

  if(verbose) {
    s = "Estimating inbreeding coefficients for %d individuals using %d markers (method: %s)"
    message(sprintf(s, length(ids), nMarkers(x), method))
  }

  if(method == "mle") {
    res = lapply(alleleData, .fest_mle)
  }
  else if(method == "simple") {
    db = getFreqDatabase(x)
    res = lapply(alleleData, .fest_simple, db = db)
  }
  else if(method == "ritland") {
    nals = nAlleles(x)
    res = lapply(alleleData, .fest_ritland, nals = nals)
  }

  res = unlist(res)

  if(verbose)
    message("Total time: ", ftime(st))

  res
}


#' @importFrom stats optimise
.fest_mle = function(dat) {

  dat = .removeMiss1(dat)

  loglik = function(f) {
    hom = dat$a1 == dat$a2
    p1 = dat$f1
    p2 = dat$f2
    n = length(hom)
    lik = ifelse(hom, f * p1 + (1-f) * p1^2, (1-f) * 2 * p1 * p2)
    sum(log(lik), na.rm = TRUE)
  }
  opt = optimise(loglik, interval = c(0,1), maximum = TRUE)
  opt$maximum
}

.fest_simple = function(dat, db) {
  if(length(dat$miss)) {
    db = db[-dat$miss]
    dat = .removeMiss1(dat)
  }
  L = length(db)
  if(L != length(dat$a1))
    stop2("Length mismatch in .fest_simple")

  # Expected homozygosity under HWE
  Hstar = sum(unlist(db)^2)

  # Observed number of homozygotes
  nHom = sum(dat$a1 == dat$a2)

  # Inbreeding estimate
  (nHom - Hstar) / (L - Hstar)
}

.fest_ritland = function(dat, nals) {
  if(length(dat$miss)) {
    nals = nals[-dat$miss]
    dat = .removeMiss1(dat)
  }
  L = length(nals)
  if(L != length(dat$a1))
    stop2("Length mismatch in .fest_ritland")

  ishom = dat$a1 == dat$a2
  numer = sum(1/dat$f1[ishom]) - L
  denom = sum(nals - 1)
  numer / denom
}

.removeMiss1 = function(dat) {
  miss = dat$miss
  if(!length(miss))
    return(dat)
  dat$a1 = dat$a1[-miss]
  dat$a2 = dat$a2[-miss]
  dat$f1 = dat$f1[-miss]
  dat$f2 = dat$f2[-miss]
  dat$miss = integer(0)
  dat
}
