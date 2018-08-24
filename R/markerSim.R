#' Marker simulation
#'
#' Simulates marker genotypes conditional on the pedigree structure and known
#' genotypes.
#'
#' This implements (with various time savers) the algorithm used in SLINK of the
#' LINKAGE/FASTLINK suite. If `partialmarker` is NULL, genotypes are simulated
#' by simple gene dropping, using [simpleSim()].
#'
#' @param x a `ped` object
#' @param N a positive integer: the number of markers to be simulated
#' @param ids a vector containing ID labels of those pedigree members whose
#'   genotypes should be simulated. By default, all individuals are included.
#' @param alleles a vector containing the alleles for the marker to be
#'   simulation. If a single integer is given, this is interpreted as the number
#'   of alleles, and the actual alleles as `1:alleles`. Must be NULL if
#'   `partialmarker` is non-NULL.
#' @param afreq a vector of length 2 containing the population frequencies for
#'   the marker alleles. Must be NULL if `partialmarker` is non-NULL.
#' @param partialmarker Either NULL (resulting in unconditional simulation), a
#'   marker object (on which the simulation should be conditioned) or the name
#'   (or index) of a marker attached to `x`.
#' @param loop_breakers a numeric containing IDs of individuals to be used as
#'   loop breakers. Relevant only if the pedigree has loops, and only if
#'   `partialmarker` is non-NULL. See [pedtools::breakLoops()].
#' @param eliminate A non-negative integer, indicating the number of iterations
#'   in the internal genotype-compatibility algorithm. Positive values can save
#'   time if `partialmarker` is non-NULL and the number of alleles is large.
#' @param seed NULL, or a numeric seed for the random number generator.
#' @param verbose a logical.
#' @return a `ped` object equal to `x` except its `markerdata` entry, which
#'   consists of the `N` simulated markers.
#' @author Magnus Dehli Vigeland
#' @seealso [simpleSim()], [pedtools::ped()]
#' @references G. M. Lathrop, J.-M. Lalouel, C. Julier, and J. Ott, *Strategies
#'   for Multilocus Analysis in Humans*, PNAS 81(1984), pp. 3443-3446.
#'
#' @examples
#' library(pedtools)
#' x = nuclearPed(2)
#' partial = marker(x, '3' = 1, alleles = 1:3)
#' markerSim(x, N = 1, alleles = 1:3)
#' markerSim(x, N = 1, partialmarker = partial)
#' markerSim(x, N = 1, partialmarker = partial)
#' markerSim(x, N = 1, ids = 4, partialmarker = partial)
#'
#' @export
markerSim = function(x, N = 1, ids = NULL, alleles = NULL, afreq = NULL, partialmarker = NULL,
  loop_breakers = NULL, eliminate = 0, seed = NULL, verbose = TRUE) {

  if (!is.ped(x) && !is.pedList(x))
      stop("x must be either a 'ped' object, a 'singleton' object, or a list of such")

  if (!is.null(seed))
      set.seed(seed)

  # if input is a list of ped objects: Apply markerSim recursively
  if (is.pedList(x))
      return(lapply(x, function(xi) markerSim(xi, N = N, ids = intersect(labels(xi), ids),
          alleles = alleles, afreq = afreq, partialmarker = partialmarker,
          loop_breakers = loop_breakers, eliminate = eliminate, verbose = verbose)))

  starttime = proc.time()

  # Reorder if necessary
  if(!has_parents_before_children(x)) {
    if(verbose) cat("Note: Changing the internal order so that all parents precede their children.\n\n")
    x = parents_before_children(x)
  }

  likel_counter = 0

  if (!is.null(x$LOOP_BREAKERS))
    stop("`ped` objects with pre-broken loops are not allowed as input to `markerSim`")
  if(is.null(ids)) ids = labels(x)

  ### Partial marker
  m = partialmarker
  if(!is.null(m)) {

    if (!is.null(alleles) || !is.null(afreq))
      stop2("When `partialmarker` is given, both 'alleles' and 'afreq' must be NULL.")

    if(is.marker(m)) {
      pedtools:::checkConsistency(x, list(m)) #TODO export this from pedtools
    }
    else if (is.atomic(m) && length(m) == 1) {
      m = getMarkers(x, markers = m)[[1]]
    }
    else
      stop2("Argument `partialmarker` must be a `marker` object, or the name (or index) of a single marker attached to `x`")

    if (is.null(attr(m, "mutmat"))) {
      err = mendelianCheck(setMarkers(x, m), verbose = F)
      if (length(err) > 0)
        stop2("Mendelian error in the given partial marker.")
    }
  }
  else {
    if (is.null(alleles))
      stop2("Marker alleles must be specified")

    if (is.numeric(alleles) && length(alleles) == 1)
      alleles = seq_len(alleles)
    m = marker(x, alleles = alleles, afreq = afreq)
  }

  alleles = alleles(m)
  afreq = unname(afreq(m))
  mutmat = attr(m, "mutmat")
  Xchrom = is_Xmarker(m)
  nall = nAlleles(m)
  mutations = allowsMutations(m)

  if (all(m == 0)) {
    return(simpleSim(x, N, alleles = alleles, afreq = afreq,
                     ids = ids, Xchrom = Xchrom,
                     mutmat = mutmat, seed = seed, verbose = verbose))
  }

  allgenos = pedprobr::allGenotypes(nall)


  gridlist = pedprobr::geno.grid.subset(x, m, labels(x), make.grid = F)


  if (verbose) {
    locus = if(!Xchrom) 'autosomal' else 'X-linked'
    plural_s = if(N>1) "s" else ""

    print(glue::glue("
      Conditional simulation of {N} {locus} marker{plural_s}.
      Target individuals: {toString(ids)}
      Conditioning on the following data:
      "))
    print(m)
  }

  # Forced genotypes:
  forcedTF = (m[, 1] == 0 | m[, 2] == 0) & (lengths(gridlist) == 1)
  m.unforced = m  # saving a copy for use in verbose output below
  for (id in (1:pedsize(x))[forcedTF])
    m[id, ] = allgenos[gridlist[[id]], ]

  if(verbose && any(forcedTF)) {
    cat("\nForced genotypes\n================\n")
    for (id in (1:pedsize(x))[forcedTF]) {
      allelchars = alleles[m[id, ]]
      if (Xchrom) allelchars = allelchars[1]
      cat(sprintf("Individual %s: %s\n", labels(x)[id], paste(allelchars, collapse = "/")))
    }
  }

  # Copies of x and m - used to determine simulation strategy - *before* possible loop breaking)
  xorig = setMarkers(x, NULL)
  morig = m

  if (loops <- x$UNBROKEN_LOOPS) {
    orig_ids = labels(x)
    x = breakLoops(setMarkers(x, m), loop_breakers = loop_breakers, verbose = verbose)
    m = x$markerdata[[1]]
    loop_breakers = x$LOOP_BREAKERS[, 1]
    gridlist = gridlist[sort.int(match(c(orig_ids, loop_breakers), orig_ids))]
  }

  ngrid = lengths(gridlist)
  SEX = x$SEX
  FIDX = x$FIDX
  MIDX = x$MIDX
  FOU = founders(x, internal=T)
  NONFOU = nonfounders(x, internal=T)

  ##### Determine simulation strategy #### Note: Using original x and m in this section (i.e.
  ##### before loop breaking)

  # Individuals that are typed (or forced - see above). Simulations condition on these.
  typedTF = (morig[, 1] != 0 | morig[, 2] != 0)
  typed = labels(xorig)[typedTF]

  # Target individuals: untyped individuals that we seek to simulate
  targets = .mysetdiff(ids, typed)
  untyped_breakers = if (loops) .mysetdiff(loop_breakers, typed) else NULL

  # Method 2: Compute joint dist of some target individuals, brute force on the remaining
  hardsim.method2 = unique.default(c(untyped_breakers, targets))
  hardsim.method2_int = internalID(x, hardsim.method2)
  method2 = .optimal.precomputation(hardsim.method2_int, N, gridlist, Xchrom, SEX = SEX)

  ### Method 3: Extend target to ancestors of targets.
  targets.plus.ancestors = .mysetdiff(c(targets, ancestors(xorig, id = targets)), typed)

  # Only ancestors of typed individuals are hard; the others can be simple-dropped
  ancestors.of.typed = ancestors(xorig, id = typed)

  hardsim.method3 = intersect(targets.plus.ancestors, ancestors.of.typed)
  hardsim.method3 = unique.default(c(untyped_breakers, hardsim.method3))

  hardsim.method3_int = internalID(x, hardsim.method3)
  method3 = .optimal.precomputation(hardsim.method3_int, N, gridlist, Xchrom, SEX = SEX)

  if (method2$calls <= method3$calls) {
    joint_int = method2$id_int
    bruteforce_int = .mysetdiff(hardsim.method2_int, joint_int)
    simpledrop = numeric()
  } else {
    joint_int = method3$id_int
    bruteforce_int = .mysetdiff(hardsim.method3_int, joint_int)
    simpledrop = .mysetdiff(targets.plus.ancestors, ancestors.of.typed)
  }

  simpledrop_int = internalID(x, simpledrop)
  simple.founders_int = intersect(simpledrop_int, FOU)
  simple.nonfounders_int = intersect(simpledrop_int, NONFOU)

  if (length(simple.nonfounders_int) > 0) {
    # Ensure sensible ordering of nonfounders (for gene dropping)
    done = c(internalID(x, typed), joint_int, bruteforce_int, simple.founders_int)
    v = simple.nonfounders_int
    v.ordered = numeric()
    while (length(v) > 0) {
      i = match(TRUE, (FIDX[v] %in% done) & (MIDX[v] %in% done))
      if (is.na(i))
        stop2("Could not determine sensible order for gene dropping. Bug report to magnusdv@medisin.uio.no is appreciated!")
      v.ordered = c(v.ordered, v[i])
      done = c(done, v[i])
      v = v[-i]
    }
    simple.nonfounders_int = v.ordered
  }

  if (verbose) {
    .printLabels = function(v) if (length(v) > 0) toString(labels(x)[v]) else "None"

    print(glue::glue("\n
      Simulation strategy
      ===================
      Pre-computed joint distribution: {.printLabels(joint_int)}
      Brute force conditional simulation: {.printLabels(bruteforce_int)}
      Hardy-Weinberg sampling (founders): {.printLabels(simple.founders_int)}
      Simple gene dropping: {.printLabels(simple.nonfounders_int)}
      Required likelihood computations:  {min(method2$calls, method3$calls)}
      \n"))
  }

  # create initial marker matrix: two columns per marker
  markers = rep.int(m, N)
  dim(markers) = c(pedsize(x), 2 * N)
  odd = seq_len(N) * 2 - 1

  if (length(joint_int) > 0) {
    allgenos_row_grid = t.default(pedprobr:::fast.grid(gridlist[joint_int]))  #Cartesian product. Each row contains 'init' row numbers of allgenos.
    jointp = apply(allgenos_row_grid, 2, function(rownrs) {
      partial = m
      partial[joint_int, ] = allgenos[rownrs, ]
      pedprobr::likelihood(x, marker1 = partial, eliminate = eliminate)
    })
    likel_counter = likel_counter + length(jointp)
    if (identical(sum(jointp), 0))
      stop("When trying to pre-compute joint probabilities: All probabilities zero. Mendelian error?")

    # fill the rows of the 'joint' individuals
    sample_rows = allgenos_row_grid[, suppressWarnings(sample.int(length(jointp), size = N,
      replace = TRUE, prob = jointp))]
    markers[joint_int, odd] = allgenos[sample_rows, 1]
    markers[joint_int, odd + 1] = allgenos[sample_rows, 2]
  }

  if (length(bruteforce_int) > 0) {
    for (i in bruteforce_int) {
      gridi = gridlist[[i]]
      rowsample = unlist(lapply(2 * seq_len(N), function(mi) {
        partial = m
        partial[] = markers[, c(mi - 1, mi)]  # preserves all attributes of the m.
        probs = unlist(lapply(gridi, function(r) {
          partial[i, ] = allgenos[r, ]
          likelihood(x, marker1 = partial, eliminate = eliminate)
        }))
        if (sum(probs) == 0) {
          print(partial)
          stop2("\nIndividual ", labels(x)[i], ": All genotype probabilities zero. Mendelian error?")
        }
        sample(gridi, size = 1, prob = probs)
      }))
      markers[i, odd] = allgenos[rowsample, 1]
      markers[i, odd + 1] = allgenos[rowsample, 2]
    }
    likel_counter = likel_counter + N * sum(ngrid[bruteforce_int])
  }

  if (length(simpledrop) > 0) {
    loopbr_int = internalID(x, x$LOOP_BREAKERS[, 1])  #integer(0) if no loops
    loopbr_dup_int = internalID(x, x$LOOP_BRREAKERS[, 2])

    if (!Xchrom)
      markers[simple.founders_int, ] = sample.int(nall, size = 2 * N * length(simple.founders_int),
        replace = TRUE, prob = afreq)
    else {
      for (f in simple.founders_int)
        markers[f, ] = switch(SEX[f],
          rep(sample.int(nall, size = N, replace = TRUE, prob = afreq), each = 2),
          sample.int(nall, size = 2 * N, replace = TRUE, prob = afreq))
    }

    markers[loopbr_dup_int, ] = markers[loopbr_int, ]  # Genotypes of the duplicated individuals. Some of these may be ungenotyped...save time by excluding these?

    for (id in simple.nonfounders_int) {
      if (!Xchrom) {
        paternal = markers[FIDX[id], odd + .rand01(N)]
        maternal = markers[MIDX[id], odd + .rand01(N)]
        if (mutations) {
          paternal = .mutate(paternal, mutmat$male)
          maternal = .mutate(maternal, mutmat$female)
        }
      } else {
        maternal = markers[MIDX[id], odd + .rand01(N)]
        if (mutations)
          maternal = .mutate(maternal, mutmat$female)

        if (SEX[id] == 1)
          paternal = maternal  # if boy, only maternal
        else {
          paternal = markers[FIDX[id], odd]  # if girl, fathers allele is forced
          if (mutations)
            paternal = .mutate(paternal, mutmat$male)
        }
      }
      markers[id, odd] = paternal
      markers[id, odd + 1] = maternal
    }
  }

  if (loops) {
    markers = markers[-x$LOOP_BREAKERS[, 2], ]
    x = tieLoops(x)
  }

  # removing genotypes for individuals that are i) originally untyped and ii) unavailable
  typedTF[forcedTF] = FALSE

  unavailable = !labels(x) %in% ids
  markers[!typedTF & unavailable, ] = 0
  attrib = attributes(partialmarker)
  attrib$name = NA_character_
  markerdata_list = lapply(seq_len(N), function(k) {
    mk = markers[, c(2 * k - 1, 2 * k)]
    attributes(mk) = attrib
    mk
  })
  class(markerdata_list) = "markerdata"
  x = setMarkers(x, markerdata_list)
  if (verbose) {
    seconds = (proc.time() - starttime)[["elapsed"]]
    print(glue::glue("
      Simulation finished
      ===================
      Number of calls to the likelihood function: {likel_counter}
      Total time used: {round(seconds, 2)} seconds
      \n"))
  }
  x
}

.mutate = function(allele_vec, mutmatrix) {
  nall = ncol(mutmatrix)
  vapply(allele_vec,
       function(a) sample.int(nall, size = 1, prob = mutmatrix[a,]),
       FUN.VALUE = 1L)
}

.optimal.precomputation = function(target_int, Nsim, gridlist, Xchrom, SEX = NULL) {
  if (length(target_int) == 0)
    return(list(calls = 0, id_int = target_int))
  if (!Xchrom) {
    nT = length(target_int)
    ngrid_target = lengths(gridlist[target_int])
    callsCum = sapply(1:nT, function(ci)
      prod(ngrid_target[seq_len(ci)]) + Nsim * sum(ngrid_target[seq_len(nT - ci) + ci]))
    minimum_index = which.min(callsCum)
    opt = list(calls = callsCum[minimum_index], id_int = target_int[seq_len(minimum_index)])
  }
  else {
    males = target_int[SEX[target_int] == 1]
    females = target_int[SEX[target_int] == 2]
    ngrid_m = lengths(gridlist[males], use.names = F)
    ngrid_f = lengths(gridlist[females], use.names = F)
    nM = length(males)
    nF = length(females)

    # Find optimal 'init' values for males/females (fewest likelihood calls)
    callsCum = matrix(nrow = nM + 1, ncol = nF + 1)
    for (ma in 0:nM) for (fe in 0:nF)
      callsCum[ma + 1, fe + 1] = prod(ngrid_m[seq_len(ma)]) * prod(ngrid_f[seq_len(fe)]) +
        Nsim * sum(c(ngrid_m[seq_len(nM - ma) + ma], ngrid_f[seq_len(nF - fe) + fe]))

    minimum_index = arrayInd(which.min(callsCum), dim(callsCum))
    id_int = c(males[seq_len(minimum_index[1] - 1)], females[seq_len(minimum_index[2] -
        1)])
    opt = list(calls = callsCum[minimum_index], id_int = id_int)
  }
  return(opt)
}






#' Unconditional marker simulation
#'
#' Unconditional simulation of unlinked markers
#'
#' This simulation is done by distributing alleles randomly to all founders,
#' followed by unconditional gene dropping down throughout the pedigree (i.e.
#' for each non-founder a random allele is selected from each of the parents).
#' Finally the genotypes of any individuals not included in `ids` are removed.
#'
#' @param x a `ped` object
#' @param N a positive integer: the number of markers to be simulated
#' @param alleles a vector containing the allele names. If missing, the alleles
#'   are taken to be `seq_along(afreq)`.
#' @param afreq a vector of length 2 containing the population frequencies for
#'   the alleles. If missing, the alleles are assumed equifrequent.
#' @param ids a vector containing ID labels of those pedigree members whose
#'   genotypes should be simulated.
#' @param Xchrom a logical: X linked markers or not?
#' @param mutmat a mutation matrix, or a list of two such matrices named
#'   'female' and 'male'. The matrix/matrices must be square, with the allele
#'   labels as dimnames, and each row must sum to 1 (after rounding to 3
#'   decimals).
#' @param seed NULL, or a numeric seed for the random number generator.
#' @param verbose a logical.
#' @return a `ped` object equal to `x` in all respects except its `markerdata`
#'   entry, which consists of the `N` simulated markers.
#' @author Magnus Dehli Vigeland
#' @seealso [markerSim()]
#'
#' @examples
#' library(pedtools)
#' x = nuclearPed(1)
#' simpleSim(x, N=3, afreq=c(0.5, 0.5))
#'
#' y = cousinPed(1, child=TRUE)
#' simpleSim(y, N=3, alleles=LETTERS[1:10])
#'
#' @importFrom utils head
#' @export
simpleSim = function(x, N, alleles, afreq, ids, Xchrom = FALSE,
                   mutmat = NULL, seed = NULL, verbose = TRUE) {
  starttime = proc.time()

  if (missing(alleles)) {
    if (missing(afreq))
      stop("Both 'alleles' and 'afreq' cannot be missing")
    alleles = seq_along(afreq)
  }

  nall = length(alleles)
  if (missing(afreq))
    afreq = rep(1, nall)/nall

  variableSNPfreqs = nall == 2 && length(afreq) != 2 && !Xchrom
  if (variableSNPfreqs)
    afreq = rep(afreq, length = N)

  if (missing(ids))
    ids = labels(x)

  mutations = !is.null(mutmat)
  if (mutations) {
    stopifnot(is.list(mutmat) || is.matrix(mutmat))
    # If single matrix given: make sex specific list
    if (is.matrix(mutmat)) {
      mutmat = pedtools:::.checkMutationMatrix(mutmat, alleles) # TODO export from pedtools
      mutmat = list(male = mutmat, female = mutmat)
    }
    else {
      stopifnot(length(mutmat) == 2, setequal(names(mutmat), c("female", "male")))
      mutmat$female = pedtools:::.checkMutationMatrix(mutmat$female, alleles)
      mutmat$male = pedtools:::.checkMutationMatrix(mutmat$male, alleles)
    }
  }

  # Reorder if necessary
  if(!has_parents_before_children(x)) {
    if(verbose) cat("Note: Changing the internal order so that all parents precede their children.\n\n")
    x = parents_before_children(x)
  }

  if (verbose) {
    locus = if(!Xchrom) 'autosomal' else 'X-linked'
    plural_s = if(N>1) "s" else ""

    print(glue::glue("
    Unconditional simulation of {N} {locus} marker{plural_s}.
    Individuals: {toString(ids)}
    "))
    if (variableSNPfreqs) {
      cat("Alleles:", toString(alleles), "\n")
      cat("Variable frequencies: p =", toString(head(afreq, 5)), ifelse(N > 5, "...\n", "\n"))
    }
    else {
      cat("Allele frequencies:\n")
      # hack to get 1 sp indent
      afr = afreq; names(afr) = alleles
      print(data.frame(as.list(afr), check.names = F), row.names = F)
    }
    cat("Mutation model:", if(mutations) "Yes" else "No", "\n\n")
  }

  if (Xchrom)
    m = .genedrop_X(x, N, nall, afreq, mutmat, seed)
  else
    m = .genedrop_AUT(x, N, nall, afreq, mutmat, seed)

  odd = seq_len(N) * 2 - 1

  m[!labels(x) %in% ids, ] = 0L
  if (variableSNPfreqs) {
    attrib = attributes(marker(x, alleles = alleles, afreq = NULL,
                               chrom = NA, mutmat = mutmat))
    frqs = as.vector(rbind(afreq, 1 - afreq))
    markerdata_list = lapply(odd, function(k) {
      mk = m[, c(k, k + 1)]
      atr = attrib
      atr$afreq = frqs[c(k, k + 1)]
      attributes(mk) = atr
      mk
    })
  } else {
    attrib = attributes(marker(x, alleles = alleles, afreq = afreq,
                               chrom = ifelse(Xchrom, 23, NA), mutmat = mutmat))
    markerdata_list = lapply(odd, function(k) {
        mk = m[, c(k, k + 1)]
        attributes(mk) = attrib
        mk
    })
  }
  x = setMarkers(x, structure(markerdata_list, class = "markerdata"))
  if (verbose) {
    seconds = (proc.time() - starttime)[["elapsed"]]
    print(glue::glue("
      Simulation finished.
      Number of calls to the likelihood function: 0.
      Total time used: {round(seconds, 2)} seconds.
      "))
  }
  x
}


.genedrop_AUT = function(x, N, nall, afreq, mutmat, seed) {
  FIDX = x$FIDX
  MIDX = x$MIDX
  SEX = x$SEX
  FOU = founders(x, internal=T)
  NONFOU = nonfounders(x, internal=T)
  mutations = !is.null(mutmat)

  m = matrix(0L, ncol = 2 * N, nrow = pedsize(x))
  odd = seq_len(N) * 2 - 1

  if (!is.null(seed)) set.seed(seed)

  variableSNPfreqs = nall==2 && length(afreq) !=2
  if (variableSNPfreqs)
    fou_alleles = unlist(lapply(afreq, function(f)
      sample.int(2, length(FOU)*2, replace = TRUE, prob = c(f, 1 - f))))
  else
    fou_alleles = sample.int(nall, size = N*length(FOU)*2, replace = TRUE, prob = afreq)

  m[FOU, ] = fou_alleles
  for (id in NONFOU) {
    paternal = m[FIDX[id], odd + .rand01(N)]
    maternal = m[MIDX[id], odd + .rand01(N)]
    if (mutations) {
      paternal = .mutate(paternal, mutmat$male)
      maternal = .mutate(maternal, mutmat$female)
    }
    m[id, odd] = paternal
    m[id, odd + 1] = maternal
  }
  m
}

.genedrop_X = function(x, N, nall, afreq, mutmat, seed) {
  FIDX = x$FIDX
  MIDX = x$MIDX
  SEX = x$SEX
  FOU = founders(x, internal=T)
  NONFOU = nonfounders(x, internal=T)
  mutations = !is.null(mutmat)

  m = matrix(0L, ncol = 2 * N, nrow = pedsize(x))
  odd = seq_len(N) * 2 - 1

  if (!is.null(seed)) set.seed(seed)

  for (f in FOU) {
    if (SEX[f] == 1)
      m[f, ] = rep(sample.int(nall, size = N, replace = TRUE, prob = afreq), each = 2)
    if (SEX[f] == 2)
      m[f, ] = sample.int(nall, size = 2 * N, replace = TRUE, prob = afreq)
  }
  for (id in NONFOU) {
    maternal = m[MIDX[id], odd + .rand01(N)]
    if (mutations)
      maternal = .mutate(maternal, mutmat$female)

    if (SEX[id] == 1) {
      paternal = maternal  # if boy, only maternal
    }
    else {
      paternal = m[FIDX[id], odd]  # if girl, fathers allele is forced
      if (mutations)
        paternal = .mutate(paternal, mutmat$male)
    }
    m[id, odd] = paternal
    m[id, odd + 1] = maternal
  }
  m
}
