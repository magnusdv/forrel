#' Read Familias .fam files
#'
#' This function parses the content of a Familias-formatted ".fam" file, and
#' converts it into suitable `ped` objects. This function does not depend on the
#' `Familias` package.
#'
#' @param famfile Path to a ".fam" file
#'
#' @return If the .fam file only contains a database, the output is a list of
#'   information (name, alleles, frequencies) about each locus. This list can be
#'   used as `locusAttributes` in e.g. [setMarkers()].
#'
#'   If the .fam file describes pedigree data, the output is a `ped` object or a
#'   list of such.
#'
#' @export
readFam = function(famfile) {
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
  message("Familias version: ", version)


  ### Individuals and genotypes

  # Number of individuals
  nid.line = if(x[4] != "") 4 else 5
  nid = as.integer(x[nid.line]) # Number of persons involved in pedigrees (but excluding "extras")
  checkInt(nid, nid.line, "number of individuals")

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
    stop(sprintf('Expected line %d to be "Known relations", but found: "%s"', id.line, x[id.line]))

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

      # Convert to data frame in ped format and insert in list
      pedigrees[[ped.idx]] = asFamiliasPedigree(id.i, fidx.i, midx.i, sex.i)
      names(pedigrees)[ped.idx] = ped.name

      # Goto next ped
      ped.line = rel.line
    }
  }

  has.probs = x[ped.line] == "#TRUE#"
  if(has.probs)
    stop("This file includes precomputed probabilities; this is not supported yet.")

  db.line = ped.line + 1
  nLoc = as.integer(x[db.line])
  checkInt(nLoc, "number of loci")

  has.info = x[db.line + 1] == "#TRUE#"
  if(has.info)
    message("Database info: ", x[db.line + 2])

  loci = vector("list", nLoc)
  loc.line = db.line + 2 + has.info
  mutWarn = T

  for(i in seq_len(nLoc)) {
    loc.name = x[loc.line]
    mutrate.fem = as.numeric(x[loc.line + 1])
    mutrate.mal = as.numeric(x[loc.line + 2])
    model.idx.fem = as.integer(x[loc.line + 3]) + 1
    model.idx.mal = as.integer(x[loc.line + 4]) + 1

    nAll.with.silent = as.integer(x[loc.line + 5]) # includes silent allele

    range.fem = as.numeric(x[loc.line + 6])
    range.mal = as.numeric(x[loc.line + 7])

    mutrate2.fem = as.numeric(x[loc.line + 8])
    mutrate2.mal = as.numeric(x[loc.line + 9])

    has.silent = x[loc.line + 10] == "#TRUE#"
    if(has.silent)
      stop("Locus ", loc.name, " has silent frequencies: this is not implemented yet")
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

    models = c("equal", "prop", "stepwise", "custom")
    if(any(c(mutrate.fem, mutrate.mal) > 0) && mutWarn) {
      warning("Mutation models are not supported yet; ignoring these", call. = F)
      mutWarn = F
    }

    # Collect locus info
    loci[[i]] = list(locusname = loc.name, alleles = frqs,
                     femaleMutationType = models[model.idx.fem],
                     femaleMutationMatrix = NULL,
                     maleMutationType = models[model.idx.mal],
                     maleMutationMatrix = NULL)

    # Goto next locus
    loc.line = loc.line + 13 + 2*nAll
  }

  ### datamatrix ###
  has.data = nid > 0 && any(lengths(sapply(datalist, '[[', "mark.idx")))
  if(!has.data)
    datamatrix = NULL
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

  if(is.null(pedigrees))
    readFamiliasLoci(loci = loci)
  else
    Familias2ped(familiasped = pedigrees, datamatrix = datamatrix, loci = loci)

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
