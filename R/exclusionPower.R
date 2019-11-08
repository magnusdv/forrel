#' Power of exclusion
#'
#' Computes the power (of a single marker) of excluding a claimed relationship,
#' given the true relationship.
#'
#' This function computes the 'Power of exclusion', as defined and discussed in
#' (Egeland et al., 2014).
#'
#' @param ped_claim a `ped` object, or a list of several such objects,
#'   describing the claimed relationship. If a list, the sets of ID labels must
#'   be disjoint, that is, all ID labels must be unique.
#' @param ped_true a `ped` object, or a list of several such objects, describing
#'   the true relationship. ID labels must be consistent with `ped_claim`.
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
#' @param verbose a logical.
#'
#' @return A single numeric value. If `plot = "plot_only"`, the function returns
#'   NULL after producing the plot.
#'
#' @author Magnus Dehli Vigeland
#'
#' @references T. Egeland, N. Pinto and M.D. Vigeland, *A general approach to
#'   power calculation for relationship testing.* Forensic Science
#'   International: Genetics 9 (2014): 186-190. doi:10.1016/j.fsigen.2013.05.001
#'
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
#' PE2 = exclusionPower(claim, true, ids = 1, alleles = als, afreq = afreq,
#'                      known_genotypes = list(c(3, 1, 1)), plot = FALSE)
#'
#' # Exclusion power if the SNP is on the X chromosome
#' PE3 = exclusionPower(claim, true, ids = ids, alleles = als, afreq = afreq,
#'                      Xchrom = TRUE, plot = FALSE)
#'
#' stopifnot(PE1 == 0.125, PE2 == 0.25, PE3 == 0.25)
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
#'                      alleles = 4, afreq = c(0.7, 0.1, 0.1, 0.1))
#'
#' stopifnot(round(c(PE1,PE2,PE3,PE4), 5) == c(0.03125, 0.00405, 0.08203, 0.03090))
#'
#' ####### How does the power change if the true pedigree is inbred?
#' sisters_LOOP = addParents(sisters, 101, father = 201, mother = 202)
#' sisters_LOOP = addParents(sisters_LOOP, 102, father = 201, mother = 203)
#'
#' # Equifrequent SNP:
#' PE5 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP,
#'                      ids = c(2, 3), alleles = 2, plot = FALSE)
#'
#' # SNP with MAF = 0.1:
#' PE6 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP,
#'                      ids = c(2, 3), alleles = 2, afreq = c(0.9, 0.1), plot = FALSE)
#'
#' stopifnot(round(c(PE5,PE6), 5) == c(0.03125, 0.00765))
#'
#' \donttest{
#' # Equifrequent tetra-allelic marker:
#' PE7 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP,
#'                      ids = c(2, 3), alleles = 4, plot = FALSE)
#'
#' # Tetra-allelic marker with one major allele:
#' PE8 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, plot = FALSE,
#'                      ids = c(2, 3), alleles = 4, afreq = c(0.7, 0.1, 0.1, 0.1))
#'
#' stopifnot(round(c(PE7,PE8), 5) == c(0.07617, 0.03457))
#' }
#'
#' @importFrom graphics abline grconvertX grconvertY layout mtext par plot rect
#' @export
exclusionPower = function(ped_claim, ped_true, ids, markerindex = NULL,
                          alleles = NULL, afreq = NULL, known_genotypes = list(),
                          Xchrom = FALSE, plot = TRUE, verbose = TRUE) {

  # st = proc.time()
  ids = as.character(ids)

  if (is.ped(ped_claim))
    ped_claim = list(ped_claim)
  if (is.ped(ped_true))
    ped_true = list(ped_true)

  # Extract ped component of each id, and check that all were found
  compsClaim = getComponent(ped_claim, ids, checkUnique = TRUE)
  if(anyNA(compsClaim))
    stop2("Individuals not found in `ped_claim`: ", ids[is.na(compsClaim)])

  # Check that all `ids` are in `true`
  compsTrue = getComponent(ped_true, ids, checkUnique = TRUE)
  if(anyNA(compsTrue))
    stop2("Individuals not found in `ped_true`: ", ids[is.na(compsTrue)])

  names(compsClaim) = names(compsTrue) = ids

  if (!is.null(markerindex)) {
    if(length(markerindex) != 1)
      stop2("Argument `markerindex` must have length 1: ", markerindex)

    ped_claim = lapply(ped_claim, function(x) selectMarkers(x, markerindex))
  }
  else {
    if (length(alleles) == 1)
      alleles = 1:alleles

    # Create and attach locus to each `claim` component
    locus = list(alleles = alleles, afreq = afreq, chrom = if (Xchrom) 23 else NA)
    ped_claim = lapply(ped_claim, function(x) setMarkers(x, locusAttributes = locus))

    # Set known genotypes
    for (tup in known_genotypes) {
      id = as.character(tup[1])
      cmp = getComponent(ped_claim, id)
      genotype(ped_claim[[cmp]], marker = 1, id = id) = tup[2:3]
    }
  }

  # Check for already typed members. TODO: Suppert for partially typed members
  typed = unlist(lapply(ped_claim, typedMembers))
  if(length(bad <- intersect(ids, typed)))
    stop2("Individual is already genotyped: ", toString(bad))

  # Transfer marker data to `true` peds
  ped_true = transferMarkers(ped_claim, ped_true)

  # Plot
  if (isTRUE(plot) || plot == "plot_only") {
    plotPedList(list(ped_claim, ped_true),
                newdev = TRUE,
                frametitles = c("Claim", "True"),
                shaded = ids,
                marker = 1)

    if (plot == "plot_only")
      return()
  }

  ###############################
  ### Computations start here ###
  ###############################

  # Step 1: Genotype combinations incompatible with "claim"
  # Step 2: Probability of incomp under `true` pedigree

  ### Step 1
  # Incompatible combinations in each component
  I.g.list = lapply(seq_along(ped_claim), function(i) {
    if(!any(compsClaim == i))
      return(NULL)

    omd = oneMarkerDistribution(ped_claim[[i]], ids = ids[compsClaim == i],
                                partialmarker = 1, verbose = FALSE, eliminate = 0)
    omd == 0
    })

  # Remove NULLs
  I.g.list = I.g.list[lengths(I.g.list) > 0]

  # outer product OR
  outerOR = function(x,y) outer(x, y, FUN = `|`)
  I.g = Reduce(outerOR, I.g.list)

  # If no incompatibilities, return 0
  if(!any(I.g))
    return(0)

  # Extract the TRUE positions (= incompatible combinations)
  incomp = which(I.g, arr.ind = TRUE, useNames = FALSE)
  colnames(incomp) = ids

  if(verbose) {
    cat("Impossible genotype combinations of `ids` given the claimed pedigree:\n")
    inc = sapply(seq_along(ids), function(i) dimnames(I.g)[[i]][incomp[, i]])
    inc = as.data.frame(inc)
    names(inc) = ids
    print(inc)
  }


  ### Step 2: Probability of incomp under `true` pedigree

  # Grid to be used as input to oneMarkerDistribution
  # Recall: Entries now refers to rows in allGenotypes(n)
  incomp.grid = incomp

  # If X: Modify male columns.  TODO: would be nice to clean this up a bit.
  m = getMarkers(ped_claim[[1]], 1)[[1]]
  if(isXmarker(m)) {
    sx = sapply(seq_along(ids), function(i) getSex(ped_true[[compsTrue[i]]], ids[i]))
    nall = length(alleles(ped_true[[1]], 1))
    if(nall > 1) {
      homoz = c(1, cumsum(nall:2) + 1)
      incomp.grid[, sx == 1] = homoz[incomp.grid[, sx == 1]]
    }
  }


  ### In the true ped: Sum probs of the incompat combinations
  p.g.list = lapply(seq_along(ped_true), function(i) {
    if(!any(compsTrue == i))
      return(NULL)

    ids.i = ids[compsTrue == i]
    grid.i = unique.matrix(incomp.grid[, ids.i, drop = FALSE])

    omd = oneMarkerDistribution(ped_true[[i]], ids.i, partialmarker = 1,
                          grid.subset = grid.i, verbose = FALSE,
                          eliminate = 1)
    omd
  })

  # Remove NULLs
  p.g.list = p.g.list[lengths(p.g.list) > 0]

  # Reduce to array
  p.g = Reduce("%o%", p.g.list)

  # The entries corresponding to exclusions
  excl = p.g[I.g]

  if(verbose) {
    cat("\nProbabilities under `true` ped:\n")
    inc = cbind(inc, prob = excl)
    print(inc)
  }

  sum(excl)
}
