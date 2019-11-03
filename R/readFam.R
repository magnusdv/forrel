#' Read Familias .fam files
#'
#' This function parses the content of a Familias-formatted ".fam" file, and
#' converts it into suitable `ped` objects. This function does not depend on the
#' `Familias` package.
#'
#' @param famfile Path to a ".fam" file.
#' @param useDVI A logical; if TRUE, the DVI section of the fam file is used to
#'   extract pedigrees and genotypes.
#' @param verbose A logical; if TRUE, various information is written to the
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
readFam = function(famfile, useDVI = F, verbose = T) {
  if(!endsWith(famfile, ".fam"))
    stop("Input file must end with '.fam'", call. = F)

  # Read entire file
  raw = readLines(famfile)
  x = gsub("\\\"", "", raw)

  # Utility function for checking integer values
  checkInt = function(a, line, txt, value) {
    if(!is.na(a)) return()
    stop(sprintf('Expected line %d to be %s, but found: "%s"',
                 line, txt, x[line]), call. = F)
  }

  # Read and print Familias version
  version = x[3]
  if(verbose)
    message("Familias version: ", version)


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

      # Add fixed relations
      nRel.i = as.integer(x[ped.line + 4])
      rel.line = ped.line + 5
      for(i in seq_len(nRel.i)) {
        par.idx = as.integer(x[rel.line]) + 1
        child.idx = as.integer(x[rel.line+1]) + 1
        if(sex[par.idx] == 1)
          fidx.i[child.idx] = par.idx
        else
          midx.i[child.idx] = par.idx

        rel.line = rel.line + 2
      }

      # Print summary
      if(verbose)
        message(sprintf(" Pedigree '%s': %d extra females, %d extra males", ped.name, nFem.i, nMal.i))

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
    dvi.families = readDVI(dvi.lines)

    if(verbose)
      message("Returning the following families:")

    famnames = names(dvi.families)
    for(i in seq_along(dvi.families)) {
      fam = dvi.families[[i]]
      if(verbose) {
        message(sprintf("%s: %d pedigrees", famnames[i], length(fam$pedigrees)))
        for(nm in names(fam$pedigrees))
          message("  ", nm)
      }
    }

    res = lapply(dvi.families, function(fam) {
      Familias2ped(familiasped = fam$pedigrees, datamatrix = fam$datamatrix,
                   loci = loci, matchLoci = T)
    })

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
  if(!is.null(pedigrees)) {
    if(verbose)
      message("\nReturning pedigrees with attached database.\n")
    Familias2ped(familiasped = pedigrees, datamatrix = datamatrix, loci = loci)
  }
  else {
    if(verbose)
      message("\nReturning database only.\n")
    readFamiliasLoci(loci = loci)
  }
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

readDVI = function(rawlines) {
  r = rawlines
  if(r[1] != "[DVI]")
    stop("I excepted the first line to be '[DVI]', but got '", r[1], "'")

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

  family_list = dvi$DVI$`Reference Families`[-1] # remove 'nFamilies'
  names(family_list) = sapply(family_list, function(fam) fam[[1]][2])

  lapply(family_list, parseFamily)
}


# Convert a "DVI Family" into a list of `datamatrix` and `pedigrees`
parseFamily = function(x) {

  persons_list = x$Persons
  if(persons_list[[c(1,1)]] == "nPersons")
    persons_list[[1]] = NULL # remove 'nPersons' entry

  ### pedigrees

  # collext id and sex of each person
  id = sapply(persons_list, function(p) p[[1]][2])
  sex = sapply(persons_list, function(p) p[[2]][2])
  sex[sex == "Male"] = 1
  sex[sex == "Female"] = 2
  sex = as.integer(sex)

  # build pedigrees defined for the family
  ped_list = x$Pedigrees[-1] # remove "nPedigrees"
  names(ped_list) = sapply(ped_list, function(pd) pd[[1]][2])

  pedigrees = lapply(ped_list, function(pd) {
    this.id = as.character(id)
    this.sex = sex

    tags = sapply(pd, '[', 1)
    vals = sapply(pd, '[', 2)

    # Data frame of parent-child pairs
    parent.tags = which(tags == "Parent")
    po = data.frame(parent = vals[parent.tags], child = vals[parent.tags + 1],
                    stringsAsFactors = F)

    # Add extra individuals if needed (e.g. "Missing person")
    if(length(extras <- setdiff(po$child, id))) {
      this.id = c(this.id, extras)
      this.sex = c(this.sex, rep(0L, length(extras)))
    }

    # Create and populate fidx and midx
    this.fidx = this.midx = integer(length(this.id))

    parent.idx = match(po$parent, this.id)
    par.is.male = this.sex[parent.idx] == 1

    child.idx = match(po$child, this.id)
    this.fidx[child.idx[par.is.male]] = parent.idx[par.is.male]
    this.midx[child.idx[!par.is.male]] = parent.idx[!par.is.male]

    # Return
    asFamiliasPedigree(this.id, this.fidx, this.midx, this.sex)
  })

  ### datamatrix
  vecs = lapply(persons_list, function(p) dnaData2vec(p$`DNA data`))
  allnames = unique(unlist(lapply(vecs, names)))
  datamatrix = do.call(rbind, lapply(vecs, function(v) v[match(allnames, names(v))]))
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
