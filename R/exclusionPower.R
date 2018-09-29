#' Power of exclusion
#'
#' Computes the power (of a single marker) of excluding a claimed relationship,
#' given the true relationship.
#'
#' This function computes the 'Power of exclusion', as defined and discussed in
#' (Egeland et al., 2014).
#'
#' @param ped_claim a `ped` object, or a list of several ped and/or singleton
#'   objects, describing the claimed relationship. If a list, the sets of ID
#'   labels must be disjoint, that is, all ID labels must be unique.
#' @param ped_true a `ped` object, or a list of several ped and/or singleton
#'   objects, describing the true relationship. ID labels must be consistent
#'   with `ped_claim`.
#' @param ids individuals available for genotyping.
#' @param markerindex NULL, or a single numeric indicating the index of a marker
#'   of `ped_claim` from which `alleles`, `afreq` and `known_genotypes` will be
#'   extracted.
#' @param alleles a numeric or character vector containing marker alleles names.
#'   Ignored if `markerindex` is non-NULL.
#' @param afreq a numerical vector with allele frequencies. An error is given if
#'   they don't sum to 1 (rounded to 3 decimals). Ignored if `markerindex` is
#'   non-NULL.
#' @param known_genotypes list of triplets `(a, b, c)`, indicating that
#'   individual `a` has genotype `b/c`. Must be NULL if `markerindex` is
#'   non-NULL.
#' @param Xchrom a logical: Is the marker on the X chromosome? Ignored if
#'   `markerindex` is non-NULL.
#' @param plot either a logical or the character 'plot_only', controlling if a
#'   plot should be produced. If 'plot_only', a plot is drawn, but no further
#'   computations are done.
#' @return A single numeric value. If `plot='plot_only'`, the function returns
#'   NULL after producing the plot.
#' @author Magnus Dehli Vigeland
#' @references T. Egeland, N. Pinto and M. D. Vigeland, *A general approach to
#'   power calculation for relationship testing.* Forensic Science
#'   International: Genetics 9 (2014): 186-190. DOI:10.1016/j.fsigen.2013.05.001
#' @examples
#'
#' ############################################
#' ### A standard case paternity case:
#' ### Compute the power of exclusion when the claimed father is in fact
#' ### unrelated to the child.
#' ############################################
#'
#' # Claim: Individual 1 is the father of indiv 3
#' claim = nuclearPed(nch = 1, sex = 2)
#'
#' # Truth: 1 and 3 are unrelated
#' true = list(singleton(id = 1), singleton(id = 3, sex = 2))
#'
#' # Indivs 1 and 3 are available for genotyping
#' ids = c(1, 3)
#'
#' # An equifrequent SNP
#' als = 2
#' afreq = c(0.5, 0.5)
#'
#' # The exclusion power assuming no known genotypes
#' PE1 = exclusionPower(claim, true, ids = ids, alleles = als, afreq = afreq)
#'
#' # Exclusion power if the child is known to have genotype 1/1:
#' PE2 = exclusionPower(claim, true, ids = ids, alleles = als, afreq = afreq,
#'                      known_genotypes = list(c(3, 1, 1)))
#'
#' # Exclusion power if the SNP is on the X chromosome
#' PE3 = exclusionPower(claim, true, ids = ids, alleles = als, afreq = afreq,
#'                      Xchrom = TRUE)
#'
#' stopifnot(PE1==0.125, PE2==0.25, PE3==0.25)
#'
#' ############################################
#' ### Example from Egeland et al. (2012):
#' ### Two females claim to be mother and daughter. Below we compute the power of various
#' ### markers to reject this claim if they in reality are sisters.
#' ############################################
#'
#' mother_daughter = nuclearPed(1, sex = 2)
#' sisters = relabel(nuclearPed(2, sex = c(2, 2)), c(101, 102, 2, 3))
#'
#' # Equifrequent SNP:
#' PE1 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, ids = c(2, 3),
#'                      alleles = 2)
#'
#' # SNP with MAF = 0.1:
#' PE2 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, ids = c(2, 3),
#'                      alleles = 2, afreq = c(0.9, 0.1))
#'
#' # Equifrequent tetra-allelic marker:
#' PE3 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, ids = c(2, 3),
#'                      alleles = 4)
#'
#' # Tetra-allelic marker with one major allele:
#' PE4 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, ids = c(2, 3),
#'                      alleles = 4, afreq=c(0.7, 0.1, 0.1, 0.1))
#'
#' stopifnot(round(c(PE1,PE2,PE3,PE4), 5) == c(0.03125, 0.00405, 0.08203, 0.03090))
#'
#' ####### How does the power change if the true pedigree is inbred?
#' sisters_LOOP = addParents(sisters, 101, father = 201, mother = 202)
#' sisters_LOOP = addParents(sisters_LOOP, 102, father = 201, mother = 203)
#'
#' # Equifrequent SNP:
#' PE5 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP,
#'                      ids = c(2, 3), alleles = 2)
#'
#' # SNP with MAF = 0.1:
#' PE6 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP,
#'                      ids = c(2, 3), alleles = 2, afreq = c(0.9, 0.1))
#'
#' stopifnot(round(c(PE5,PE6), 5) == c(0.03125, 0.00765))
#'
#' \dontrun{
#' # Equifrequent tetra-allelic marker:
#' PE7 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP,
#'                      ids = c(2, 3), alleles = 4)
#'
#' # Tetra-allelic marker with one major allele:
#' PE8 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP,
#'                      ids = c(2, 3), alleles = 4, afreq = c(0.7, 0.1, 0.1, 0.1))
#'
#' stopifnot(round(c(PE7,PE8), 5) == c(0.07617, 0.03457))
#' }
#'
#' @importFrom graphics abline grconvertX grconvertY layout mtext par plot rect
#' @export
exclusionPower = function(ped_claim, ped_true, ids, markerindex = NULL,
                          alleles = NULL, afreq = NULL, known_genotypes = list(),
                          Xchrom = FALSE, plot = TRUE) {
  if (is.ped(ped_claim))
    ped_claim = list(ped_claim)
  if (is.ped(ped_true))
    ped_true = list(ped_true)

  ids_claim = lapply(ped_claim, function(x)
    ids[ids %in% labels(x)]) #internalID?
  ids_true = lapply(ped_true, function(x)
    ids[ids %in% labels(x)])

  N_claim = length(ped_claim)
  N_true = length(ped_true)
  N = N_claim + N_true

  if (is.null(alleles)) {
    # Use markerdata of ped_claim and ped_true.  NB: No compatibility testing is done!!
    partial_claim = lapply(ped_claim, function(p)
      p$markerdata[[markerindex]])
    partial_true = lapply(ped_true, function(p)
      p$markerdata[[markerindex]])
  } else {
    if (length(alleles) == 1)
      alleles = 1:alleles

    chrom = if (Xchrom) 23 else NA

    partial_claim = lapply(1:N_claim, function(i) {
      x = ped_claim[[i]]
      m = marker(x,
                 alleles = alleles,
                 afreq = afreq,
                 chrom = chrom)
      for (tup in known_genotypes) {
        id = tup[1]
        if (id %in% labels(x))
          genotype(m, id) = tup[2:3]
      }
      m
    })
    partial_true = lapply(1:N_true, function(i) {
      x = ped_true[[i]]
      m = marker(x,
                 alleles = alleles,
                 afreq = afreq,
                 chrom = chrom)
      for (tup in known_genotypes) {
        id = tup[1]
        if (id %in% labels(x))
          genotype(m, id) = tup[2:3]
      }
      m
    })
  }

  if (isTRUE(plot) || plot == "plot_only") {
    op = par(oma = c(0, 0, 3, 0), xpd = NA)
    widths = ifelse(sapply(c(ped_claim, ped_true), is.singleton), 1, 2)
    claim_ratio = sum(widths[1:N_claim]) / sum(widths)
    layout(rbind(1:N), widths = widths)
    has_genotypes = length(known_genotypes) > 0
    for (i in 1:N) {
      if (i <= N_claim) {
        x = ped_claim[[i]]
        avail = ids_claim[[i]]
        mm = if (has_genotypes) partial_claim[[i]] else NULL
      } else {
        x = ped_true[[i - N_claim]]
        avail = ids_true[[i - N_claim]]
        mm = if (has_genotypes) partial_true[[i - N_claim]] else NULL
      }

      cols = ifelse(labels(x) %in% avail, 2, 1)
      plot(x, marker = mm, col = cols, margin = c(2, 4, 2, 4), title = "")
    }

    mtext("Claim", outer = TRUE, at = claim_ratio / 2)
    mtext("True", outer = TRUE, at = 0.5 + claim_ratio / 2)
    rect(
      grconvertX(0.02, from = "ndc"),
      grconvertY(0.02, from = "ndc"),
      grconvertX(claim_ratio - 0.02, from = "ndc"),
      grconvertY(0.98, from = "ndc")
    )
    rect(
      grconvertX(claim_ratio + 0.02, from = "ndc"),
      grconvertY(0.02, from = "ndc"),
      grconvertX(0.98, from = "ndc"),
      grconvertY(0.98, from = "ndc")
    )
    par(op)
    if (plot == "plot_only")
      return()
  }

  p.g = Reduce("%o%", lapply(which(lengths(ids_true) > 0), function(i)
    oneMarkerDistribution(
      ped_true[[i]],
      ids = ids_true[[i]],
      partialmarker = partial_true[[i]],
      verbose = F
    )))
  I.g = Reduce("%o%", lapply(which(lengths(ids_claim) > 0), function(i)
    oneMarkerDistribution(
      ped_claim[[i]],
      ids = ids_claim[[i]],
      partialmarker = partial_claim[[i]],
      verbose = F
    ) == 0))

  sum(p.g * I.g)
}

