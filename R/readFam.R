#' Read `Familias` .fam files
#'
#' This function parses the content of a `Familias`-formatted ".fam" file, and
#' converts it into suitable `ped` objects. This function does not depend on the
#' `Familias` R package.
#'
#' @param famfile Path to a ".fam" file.
#' @param useDVI A logical, indicating if the DVI section of the fam file should
#'   be identified and parsed. If `NA` (the default), the DVI section is
#'   included if it is present in the input file.
#' @param Xchrom A logical. If TRUE, the `chrom` attribute of all markers will
#'   be set to "X". (Default = FALSE.)
#' @param verbose A logical. If TRUE, various information is written to the
#'   screen during the parsing process.
#'
#' @return If the .fam file only contains a database, the output is a list of
#'   information (name, alleles, frequencies) about each locus. This list can be
#'   used as `locusAttributes` in e.g. [setMarkers()].
#'
#'   If the .fam file describes pedigree data, the output is a `ped` object or a
#'   list of such.
#'
#'   If `useDVI = TRUE`, then the families described under `Reference Families`
#'   are parsed and converted to `ped` objects. Each family generally describes
#'   multiple pedigrees, so the output gets another layer in this case.
#'
#' @importFrom pedmut mutationMatrix
#' @export
readFam = function(famfile, useDVI = NA, Xchrom = FALSE, verbose = TRUE) {
  if(!endsWith(famfile, ".fam"))
    stop("Input file must end with '.fam'", call. = FALSE)

  # Read entire file
  raw = readLines(famfile)
  x = gsub("\\\"", "", raw)

  # Utility function for checking integer values
  checkInt = function(a, line, txt, value) {
    if(!is.na(a)) return()
    stop(sprintf('Expected line %d to be %s, but found: "%s"',
                 line, txt, x[line]), call. = FALSE)
  }

  # Read and print Familias version
  version = x[3]
  if(verbose)
    message("Familias version: ", version)

  if(is.na(useDVI))
    useDVI = "[DVI]" %in% x
  else if(useDVI && !"[DVI]" %in% x)
    stop2("`useDVI` is TRUE, but no line equals '[DVI]'")

  if(verbose)
    message("Read DVI: ", if(useDVI) "Yes" else "No")

  ### Individuals and genotypes

  # Number of individuals
  nid.line = if(x[4] != "") 4 else 5
  nid = as.integer(x[nid.line]) # Number of persons involved in pedigrees (but excluding "extras")
  checkInt(nid, nid.line, "number of individuals")
  if(verbose)
    message("\nNumber of individuals (excluding 'extras'): ", nid)

  # Initialise id & sex
  id = character(nid)
  sex = integer(nid)

  # Initialise list holding genotypes (as allele indices)
  datalist = vector("list", nid)

  # Read data for each individual
  id.line = nid.line + 1
  for(i in seq_len(nid)) {
    id[i] = x[id.line]
    sex[i] = ifelse(x[id.line + 4] == "#TRUE#", 1, 2)

    nmi = as.integer(x[id.line + 5])
    checkInt(nmi, id.line + 5,
             sprintf('number of genotypes for "%s"', id[i]))
    if(verbose)
      message(sprintf("  Individual '%s': Genotypes for %d markers read", id[i], nmi))

    a1.lines = seq(id.line + 6, by = 3, length = nmi)
    a1.idx = as.integer(x[a1.lines]) + 1
    a2.idx = as.integer(x[a1.lines + 1]) + 1
    mark.idx = as.integer(x[a1.lines + 2]) + 1
    datalist[[i]] = list(id = id[i], a1.idx = a1.idx,
                         a2.idx = a2.idx, mark.idx = mark.idx)

    id.line = id.line + 6 + 3*nmi
  }

  ### Fixed relations

  # Storage for twins
  twins = list()

  kr.line = id.line
  if(x[kr.line] != "Known relations")
    stop2(sprintf('Expected line %d to be "Known relations", but found: "%s"', id.line, x[id.line]))

  # Add extras to id & sex
  nFem = as.integer(x[kr.line + 1])
  nMal = as.integer(x[kr.line + 2])
  id = c(id, character(length = nFem + nMal))
  sex = c(sex, rep.int(2:1, c(nFem, nMal)))

  # Initialise fidx, midx
  fidx = midx = integer(length(id))

  # Add fixed relations
  nRel = as.integer(x[kr.line + 3])
  rel.line = kr.line + 4
  for(i in seq_len(nRel)) {
    par.idx = as.integer(x[rel.line]) + 1
    child.idx = as.integer(x[rel.line+1]) + 1
    if(sex[par.idx] == 1)
      fidx[child.idx] = par.idx
    else
      midx[child.idx] = par.idx

    # Goto next relation
    rel.line = rel.line + 2
  }

  # Initialise list of final pedigrees
  nPed = as.integer(x[rel.line])
  checkInt(nPed, "number of pedigrees")
  if(verbose)
    message("\nNumber of pedigrees: ", nPed)

  # If no more pedigree info, finish pedigree part
  if(nPed == 0) {
    pedigrees = asFamiliasPedigree(id, fidx, midx, sex)
  }
  ped.line = rel.line + 1

  ### Additional relationships, unique to each ped
  if(nPed > 0) {
    pedigrees = vector("list", nPed)

    # Process each pedigree
    for(i in seq_len(nPed)) {
      ped.idx = as.integer(x[ped.line]) + 1
      ped.name = x[ped.line + 1]

      # Add extras in the i'th pedigree
      nFem.i = as.integer(x[ped.line + 2])
      nMal.i = as.integer(x[ped.line + 3])
      id.i = c(id, character(nFem.i + nMal.i))
      sex.i = c(sex, rep.int(2:1, c(nFem.i, nMal.i)))
      fidx.i = c(fidx, integer(nFem.i + nMal.i))
      midx.i = c(midx, integer(nFem.i + nMal.i))

      # Print summary
      if(verbose)
        message(sprintf(" Pedigree '%s' (%d extra females, %d extra males)", ped.name, nFem.i, nMal.i))

      # Add fixed relations
      nRel.i = as.integer(x[ped.line + 4])
      rel.line = ped.line + 5
      for(i in seq_len(nRel.i)) {
        par.idx = as.integer(x[rel.line]) + 1
        child.idx = as.integer(x[rel.line+1]) + 1
        if(is.na(par.idx)) {
          if(grepl("Direct", x[rel.line])) {
            par.idx = as.integer(substring(x[rel.line], 1, 1)) + 1
            twins = c(twins, list(par.idx, child.idx))
            if(verbose) message("  Twins: ", toString(id[c(par.idx, child.idx)]))
            stop2("File contains twins - this is not supported yet")
          }
        }
        if(sex[par.idx] == 1)
          fidx.i[child.idx] = par.idx
        else
          midx.i[child.idx] = par.idx

        rel.line = rel.line + 2
      }


      # Convert to familiaspedigree and insert in list
      pedigrees[[ped.idx]] = asFamiliasPedigree(id.i, fidx.i, midx.i, sex.i)
      names(pedigrees)[ped.idx] = ped.name

      # Goto next ped
      ped.line = rel.line
    }
  }

  has.probs = x[ped.line] == "#TRUE#"
  if(has.probs)
    stop("\nThis file includes precomputed probabilities; this is not supported yet.")

  ### Database ###

  db.line = ped.line + 1
  nLoc = as.integer(x[db.line])
  checkInt(nLoc, "number of loci")

  has.info = x[db.line + 1] == "#TRUE#"
  if(verbose) {
    if(has.info)
      message("\nDatabase: ", x[db.line + 2])
    else
      message("")
  }

  if(verbose)
    message("Number of loci: ", nLoc)

  loci = vector("list", nLoc)
  loc.line = db.line + 2 + has.info

  # Storage for unsupported mutation models
  unsupp = character()


  for(i in seq_len(nLoc)) {
    loc.name = x[loc.line]
    mutrate.fem = as.numeric(x[loc.line + 1])
    mutrate.mal = as.numeric(x[loc.line + 2])
    model.idx.fem = as.integer(x[loc.line + 3])
    model.idx.mal = as.integer(x[loc.line + 4])

    nAll.with.silent = as.integer(x[loc.line + 5]) # includes silent allele

    range.fem = as.numeric(x[loc.line + 6])
    range.mal = as.numeric(x[loc.line + 7])

    mutrate2.fem = as.numeric(x[loc.line + 8])
    mutrate2.mal = as.numeric(x[loc.line + 9])

    has.silent = x[loc.line + 10] == "#TRUE#"
    if(has.silent)
      stop2("Locus ", loc.name, " has silent frequencies: this is not implemented yet")
    silent.freq = as.numeric(x[loc.line + 11])

    # Number of alleles except the silent
    nAll = x[loc.line + 12]
    nAll = as.integer(strsplit(nAll, "\t")[[1]][[1]])
    checkInt(nAll, sprintf('number of alleles for marker "%s"', loc.name))

    # Read alleles and freqs
    als.lines = seq(loc.line + 13, by = 2, length.out = nAll)
    als = x[als.lines]
    frqs = as.numeric(x[als.lines + 1])
    names(frqs) = als

    # Mutation models
    models = c("equal", "proportional", "stepwise", "3", "4")
    maleMod = models[model.idx.mal + 1]
    femaleMod = models[model.idx.fem + 1]

    if(maleMod %in% models[1:3]) {
      maleMutMat = mutationMatrix(model = maleMod, alleles = als, afreq = frqs,
                                  rate = mutrate.mal, rate2 = mutrate2.mal)
    }
    else {
      if(verbose)
        message(sprintf("*** Ignoring male mutation model '%s' ***", maleMod))
      maleMutMat = NULL
      unsupp = c(unsupp, maleMod)
    }

    if(femaleMod %in% models[1:3]) {
      femaleMutMat = mutationMatrix(model = femaleMod, alleles = als, afreq = frqs,
                                    rate = mutrate.fem, rate2 = mutrate2.fem)
    }
    else {
      if(verbose)
        message(sprintf("*** Ignoring female mutation model '%s' ***", femaleMod))
      femaleMutMat = NULL
      unsupp = c(unsupp, femaleMod)
    }

    # Print locus summary
    if(maleMod == femaleMod && mutrate.fem == mutrate.mal)
      mut_txt = sprintf("model = %s, rate = %.2g (unisex)", maleMod, mutrate.mal)
    else
      mut_txt = sprintf("male model = %s, male rate = %.2g, female model = %s, female rate = %.2g",
                         maleMod, mutrate.mal, femaleMod, mutrate.fem)
    if(verbose) message(sprintf("  %s: %d alleles, %s", loc.name, length(frqs), mut_txt))

    # Collect locus info
    loci[[i]] = list(locusname = loc.name, alleles = frqs,
                     femaleMutationType = femaleMod,
                     femaleMutationMatrix = femaleMutMat,
                     maleMutationType = maleMod,
                     maleMutationMatrix = maleMutMat)

    # Goto next locus
    loc.line = loc.line + 13 + 2*nAll
  }

  # Warn about unsupported mutation models
  if(length(unsupp) > 0) {
    unsupp = sort(unique.default(unsupp))
    tmp ="Some mutation models were set to NULL. Model '%s' is not supported yet."
    warning(paste0(sprintf(tmp, unsupp), collapse = "\n"), call. = FALSE)
  }


  ###########
  ### DVI ###
  ###########

  if(useDVI) {
    if(verbose)
      message("\n*** Reading DVI section ***")
    dvi.start = match("[DVI]", raw)
    if(is.na(dvi.start))
      stop2("Expected keyword '[DVI]' not found")
    dvi.lines = raw[dvi.start:length(raw)]
    dvi.families = readDVI(dvi.lines, verbose = verbose)

    if(verbose)
      message("*** Finished DVI section ***\n")

    if(verbose)
      message("\nConverting to `ped` format")
    res = lapply(dvi.families, function(fam) {
      Familias2ped(familiasped = fam$pedigrees, datamatrix = fam$datamatrix,
                   loci = loci, matchLoci = TRUE)
    })

    # Set all chrom attributes to X if indicated
    if(Xchrom) {
      if(verbose) message("Changing all chromosome attributes to `X`")
      chrom(res, seq_along(loci)) = "X"
    }

    if(verbose) message("")
    return(res)
  }

  ##################
  ### If not DVI ###
  ##################

  ### datamatrix ###
  has.data = nid > 0 && any(lengths(sapply(datalist, '[[', "mark.idx")))
  if(!has.data) {
    datamatrix = NULL
  }
  else {
    loc.names = vapply(loci, function(ll) ll$locusname, FUN.VALUE = "")

    # Organise observed alleles in two index matrices
    dm.a1.idx = dm.a2.idx = matrix(NA, nrow = nid, ncol = nLoc,
                                   dimnames = list(id[seq_len(nid)], loc.names))
    for(i in seq_len(nid)) {
      g = datalist[[i]]
      dm.a1.idx[i, g$mark.idx] = g$a1.idx
      dm.a2.idx[i, g$mark.idx] = g$a2.idx
    }

    # Initalise data matrix
    dmn = list(id[seq_len(nid)], paste(rep(loc.names, each = 2), 1:2, sep = "."))
    datamatrix = matrix(NA_character_, nrow = nid, ncol = 2*nLoc, dimnames = dmn)

    # Fill in observed alleles
    for(i in seq_len(nLoc)) {
      als.i = names(loci[[i]]$alleles)
      datamatrix[, 2*i - 1] = als.i[dm.a1.idx[, i]]
      datamatrix[, 2*i]     = als.i[dm.a2.idx[, i]]
    }
  }

  # Return
  if(!is.null(pedigrees)) {
    if(verbose)
      message("\nConverting to `ped` format")
    res = Familias2ped(familiasped = pedigrees, datamatrix = datamatrix, loci = loci)

    # Set all chrom attributes to X if indicated
    if(Xchrom) {
      if(verbose) message("Changing all chromosome attributes to `X`")
      chrom(res, seq_along(loci)) = "X"
    }
  }
  else {
    if(verbose)
      message("\nReturning database only")
    res = readFamiliasLoci(loci = loci)
  }

  if(verbose) message("")
  res
}




# Create a FamiliasPedigree from scratch (without loading Familias package)
asFamiliasPedigree = function(id, findex, mindex, sex) {
  if(length(id) == 0)
    return(NULL)

  if(is.numeric(sex))
    sex = ifelse(sex == 1, "male", "female")

  x = list(id = id, findex = findex, mindex = mindex, sex = sex)
  class(x) = "FamiliasPedigree"

  x
}

#########################################
### Utilities for parsing DVI section ###
#########################################

readDVI = function(rawlines, verbose = TRUE) {
  r = rawlines
  if(r[1] != "[DVI]")
    stop("Expected the first line of DVI part to be '[DVI]', but got '", r[1], "'")

  ### Parse raw lines into nested list named `dvi`
  dvi = list()
  ivec = character()

  # number of brackets on each line
  brackets = as.integer(regexpr("[^[]", r)) - 1

  # pre-split lines
  splits = strsplit(r, "= ")

  # Populate `dvi` list
  for(i in seq_along(r)) {
    line = r[i]
    if(line == "")
      next
    br = brackets[i]

    if(br == 0) {
      dvi[[ivec]] = c(dvi[[ivec]], list(splits[[i]]))
    }
    else {
      name = gsub("[][]", "", line)
      ivec = c(ivec[seq_len(br - 1)], name)
      dvi[[ivec]] = list()
    }

  }

  # Initialise output list
  res = list()

  # Unidentified persons, if any
  un = parseUnidentified(dvi$DVI$`Unidentified persons`, verbose = verbose)
  if(!is.null(un))
    res$`Unidentified persons` = un

  # Reference families
  refs_raw = dvi$DVI$`Reference Families`
  refs = refs_raw[-1] # remove 'nFamilies'
  stopifnot((nFam <- length(refs)) == as.integer(refs_raw[[c(1,2)]]))
  if(verbose)
    message("\nReference families: ", nFam)

  names(refs) = sapply(refs, function(fam) fam[[1]][2])
  refs = lapply(refs, parseFamily, verbose = verbose)

  # Return
  c(res, refs)
}

parseUnidentified = function(x, verbose = TRUE) {
  if(length(x) == 0)
    return(NULL)

  nPers = x[[c(1,2)]]
  if(verbose)
    message("Unidentified persons: ", nPers)

  if(nPers == "0")
    return(NULL)

  x = x[-1]

  ### id and sex
  id = sapply(x, function(p) getValue(p[[1]], iftag = "Name", NA))
  sex = sapply(x, function(p) getValue(p[[2]], iftag = "Gender", 0))
  sex[sex == "Male"] = 1
  sex[sex == "Female"] = 2
  s = asFamiliasPedigree(as.character(id), 0, 0, as.integer(sex))

  if(verbose)
    for(nm in id) message("  ", nm)

  ### datamatrix
  vecs = lapply(x, function(p) dnaData2vec(p$`DNA data`))

  # Remove NULLs
  vecs = vecs[!sapply(vecs, is.null)]

  # All column names
  allnames = unique(unlist(lapply(vecs, names)))

  # Ensure same order in each vector, and fill in NA's
  vecs_ordered = lapply(vecs, function(v) structure(v[allnames], names = allnames))

  # Bind to matrix
  datamatrix = do.call(rbind, vecs_ordered)
  rownames(datamatrix) = id[rownames(datamatrix)]

  ### return
  list(pedigrees = s, datamatrix = datamatrix)
}

# Convert a "DVI Family" into a list of `datamatrix` and `pedigrees`
parseFamily = function(x, verbose = TRUE) {

  famname = x[[c(1,2)]]
  nPers = as.integer(x$Persons[[c(1,2)]])
  nPeds = as.integer(x$Pedigrees[[c(1,2)]])

  if(verbose)
    message(sprintf("  %s (%d persons, %d pedigrees)", famname, nPers, nPeds))

  ### Persons
  persons_list = x$Persons[-1]

  id = sapply(persons_list, function(p) getValue(p[[1]], iftag = "Name", NA))
  sex = sapply(persons_list, function(p) getValue(p[[2]], iftag = "Gender", 0))

  sex[sex == "Male"] = 1
  sex[sex == "Female"] = 2
  sex = as.integer(sex)

  ### pedigrees
  ped_list = x$Pedigrees[-1] # remove "nPedigrees"
  names(ped_list) = sapply(ped_list, function(pd) getValue(pd[[1]], iftag = "Name", NA))

  pedigrees = lapply(ped_list, function(pd) {
    if(verbose)
      message("    ", pd[[1]][2])

    this.id = as.character(id)
    this.sex = sex

    tags = sapply(pd, '[', 1)
    vals = sapply(pd, '[', 2)

    # Data frame of parent-child pairs
    parent.tags = which(tags == "Parent")
    po = data.frame(parent = vals[parent.tags], child = vals[parent.tags + 1],
                    stringsAsFactors = FALSE)

    # Add extra individuals if needed (e.g. "Missing person")
    extras <- setdiff(c(po$parent, po$child), id)
    if(length(extras)) {
      this.id = c(this.id, extras)
      this.sex = c(this.sex, rep(0L, length(extras)))
    }

    names(this.sex) = this.id
    parent.sex = this.sex[po$parent]

    # Try to fix parents with undecided sex
    if(any(parent.sex == 0)) {
      par.nosex = unique(po$parent[parent.sex == 0])
      for(p in par.nosex) {
        chi = po$child[po$parent == p] # children of him/her
        spou = unique(setdiff(po$parent[po$child %in% chi], p))
        if(all(this.sex[spou] == 1))
          this.sex[p] = 2
        else if(all(this.sex[spou] == 2))
          this.sex[p] = 1
        else
          stop2("Cannot decide sex of this parent: ", p)
      }

      # Now try again
      parent.sex = this.sex[po$parent]
    }

    parent.idx = match(po$parent, this.id)
    child.idx = match(po$child, this.id)

    # Create and populate fidx and midx
    this.fidx = this.midx = integer(length(this.id))
    this.fidx[child.idx[parent.sex == 1]] = parent.idx[parent.sex == 1]
    this.midx[child.idx[parent.sex == 2]] = parent.idx[parent.sex == 2]

    # Return
    asFamiliasPedigree(this.id, this.fidx, this.midx, this.sex)
  })

  ### datamatrix
  vecs = lapply(persons_list, function(p) dnaData2vec(p$`DNA data`))

  # Remove NULLs
  vecs = vecs[!sapply(vecs, is.null)]

  # All column names
  allnames = unique(unlist(lapply(vecs, names)))

  # Ensure same order in each vector, and fill in NA's
  vecs_ordered = lapply(vecs, function(v) structure(v[allnames], names = allnames))

  # Bind to matrix
  datamatrix = do.call(rbind, vecs_ordered)
  rownames(datamatrix) = id[rownames(datamatrix)]

  ### return
  list(pedigrees = pedigrees, datamatrix = datamatrix)
}


# DNA data for single person --> named vector
dnaData2vec = function(x) {
  dat = do.call(rbind, x)
  val = dat[, 2]

  idx = which(dat[,1] == "SystemName")
  nLoc = length(idx)
  if(nLoc == 0)
    return()

  res = character(2 * nLoc)
  res[2*(1:nLoc) - 1] = val[idx + 1]
  res[2*(1:nLoc)] = val[idx + 2]
  names(res) = paste(rep(val[idx], each = 2), 1:2, sep = ".")
  res
}

getValue = function(x, iftag, default) {
  if(x[1] == iftag) x[2] else default
}
