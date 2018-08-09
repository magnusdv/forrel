#' Convert `Familias` objects to `ped` format
#'
#' Familias is a widely used program for computations in forensic genetics. The
#' function documented here converts pedigrees and marker data from Familias to
#' [pedtools::ped()] format, thus enabling such data to be analysed with
#' `forrel`. This may be of interest for specialized computations not
#' implemented in Familias, e.g. conditional simulations.
#'
#' The Familias program represents pedigrees and marker data in a way that
#' differs from `pedtools` in several ways, mostly because of the latter's
#' stricter definition of a *pedigree*. In `pedtools` pedigrees must be
#' connected, and each member must have either 0 or 2 parents present in the
#' pedigree. None of this is required by `FamiliasPedigree` objects. The
#' conversion function `Familias2ped` takes care of all potential differences:
#' It converts each FamiliasPedigree into a list of connected ped objects,
#' adding missing parents where needed.
#'
#' @param familiasped A [Familias::FamiliasPedigree()] object or a list of such.
#' @param datamatrix A data frame with two columns per marker (one for each
#'   allele) and one row per individual.
#' @param loci A [Familias::FamiliasLocus()] object or a list of such.
#'
#' @return A `ped` object, or a list of such.
#' @author Magnus Dehli Vigeland, Thore Egeland
#'
#' @references Windows Familias is freely availabe from <http://familias.name>.
#' @examples
#'
#' \dontrun{
#' # Example 1
#' # =========
#' library(Familias)
#' library(pedtools)
#' data(NorwegianFrequencies)
#' TH01 = NorwegianFrequencies$TH01
#' locus1 = FamiliasLocus(TH01)
#' persons = c('mother', 'daughter', 'AF')
#' ped1 = FamiliasPedigree(id = persons,
#'                         dadid = c(NA, 'AF', NA),
#'                         momid = c(NA, 'mother', NA),
#'                         sex = c('female', 'female', 'male'))
#' datamatrix = data.frame(THO1.1=c(NA, 8, NA), THO1.2=c(NA,9.3, NA))
#' rownames(datamatrix) = persons
#' x = Familias2ped(ped1, datamatrix, locus1)
#' plot(x, marker = 1)
#'
#' # Example 2 TODO - test!
#' # =========
#' library(fam2r)
#' data(F21)
#' pedigrees = F21$pedigrees
#' datamatrix = F21$datamatrix
#' loci = F21$loci
#'
#' x = Familias2ped(pedigrees, datamatrix, loci)
#' plotPedList(x, new=TRUE, frametitles=c('H1', 'H2'))
#'
#' # Give dev.width explicitly to allow for long names
#' plotPedList(x, new = TRUE, frametitles = c('H1', 'H2'), dev.width = 17)
#'
#' # Numerical labels work better
#' plotPedList(x, new=TRUE, id.labels = 'num', frametitles = c('H1', 'H2'))
#'
#' }
#'
#' @export
Familias2ped = function(familiasped, datamatrix, loci) {

  ### If first argument is a list of FamiliasPedigrees, convert one at a time.
  if (is.list(familiasped) && class(familiasped[[1]]) == "FamiliasPedigree") {
      res = lapply(familiasped, function(p)
        Familias2ped(p, datamatrix = datamatrix, loci = loci))
      return(res)
  }

  ### Part 1: pedigree
  id = familiasped$id
  findex = familiasped$findex
  mindex = familiasped$mindex
  p = data.frame(id = id,
                 fid = ifelse(findex > 0, id[findex], 0),
                 mid = ifelse(mindex > 0, id[mindex], 0),
                 sex = ifelse(familiasped$sex == "male", 1, 2),
                 stringsAsFactors = FALSE)

  fatherMissing = which(p$fid == 0 & p$mid != 0)
  motherMissing = which(p$fid != 0 & p$mid == 0)

  nFath = length(fatherMissing)
  nMoth = length(motherMissing)

  newFathers = paste("added", seq(1, length.out=nFath), sep="_")
  newMothers = paste("added", seq(nFath + 1, length.out=nMoth), sep="_")

  # add new fathers
  if (nFath > 0) {
    p = rbind(p, data.frame(id = newFathers, fid = 0, mid = 0, sex = 1))
    p[fatherMissing, "fid"] = newFathers
  }

  # add new mothers
  if (nMoth > 0) {
    p = rbind(p, data.frame(id = newMothers, fid = 0, mid = 0, sex = 2))
    p[motherMissing, "mid"] = newMothers
  }

  # identify connected components and add famid column. Keep order unchanged!
  comps = connectedComponents(p$id, p$fid, p$mid)
  famid = integer(nrow(p))
  for (i in 1:length(comps)) famid[match(comps[[i]], p$id)] = i

  p = cbind(famid = famid, p)

  ### Part2: datamatrix

  if (!is.null(datamatrix)) {
    # sort datamatrix according to ped order
    id_idx = match(familiasped$id, rownames(datamatrix))
    if (anyNA(id_idx))
      stop2("ID label not found among the datamatrix rownames: ", setdiff(familiasped$id, rownames(datamatrix)))
    datamatrix = datamatrix[id_idx, , drop = FALSE]

    # convert from factor to character
    allelematrix = sapply(datamatrix, as.character)

    # replace NA with 0
    allelematrix[is.na(allelematrix)] = "0"

    # add empty rows corresponding to new parents
    addedParents = matrix("0", nrow = nFath + nMoth, ncol = ncol(allelematrix))
    allelematrix = rbind(allelematrix, addedParents)

    p = cbind(p, allelematrix, stringsAsFactors = FALSE)
  }

  ### Part 3: marker annotations

  annotations = readFamiliasLoci(loci)

  ### Create ped object

  pedtools::as.ped(p, locus_annotations = annotations)
}




#' @export
#' @rdname Familias2ped
readFamiliasLoci = function(loci) {
  if (is.null(loci))
    return(NULL)
  if (class(loci) == "FamiliasLocus")
    loci = list(loci)

  lapply(loci, function(a) {
    malemut = a$maleMutationMatrix
    femalemut = a$femaleMutationMatrix

    if (all(diag(malemut) == 1)) malemut = NULL
    if (all(diag(femalemut) == 1)) femalemut = NULL

    if (is.null(malemut) && is.null(femalemut))
      mutmat = NULL
    else
      mutmat = list(male = malemut, female = femalemut)

    list(name = a$locusname, alleles = names(a$alleles),
         afreq = as.numeric(a$alleles), mutmat = mutmat)
  })
}



#' Connected pedigree components
#'
#' Compute the connected parts of a (possibly) disconnected pedigree. This is a
#' necessary step when converting pedigree data from software like Familias
#' (which allows disconnected pedigrees) to `pedtools` (which requires pedigrees
#' to be connected).
#'
#' @param ID A vector of ID labels (character or numeric)
#' @param FID The ID labels of the fathers (or `0` if missing)
#' @param MID The ID labels of the mothers (or `0` if missing)
#'
#' @return A list, where each element is a vector of ID labels constituting a
#'   connected pedigree
#' @export
connectedComponents = function(ID, FID, MID) {
  # Placeholder for final components
  comps = list()

  # Starting point: List of all founders and trios
  temp = lapply(seq_along(ID), function(i) .mysetdiff(c(ID[i], FID[i], MID[i]), 0))  # trios

  while (length(temp) > 0) {
    # Check if first vector overlaps with any of the others
    a = temp[[1]]
    remov = numeric()
    for (j in seq_along(temp)[-1]) {
      if (any(match(a, temp[[j]], nomatch = 0) > 0)) {
        a = unique.default(c(a, temp[[j]]))
        remov = c(remov, j)
      }
    }

    if (length(remov) > 0) {
      # Remove any overlapping vectors, and insert the union as first element
      temp[remov] = NULL
      temp[[1]] = a
    } else {
      # If no overlaps, we have a maximal component. Move to comps and remove from temp.
      comps = c(comps, list(sort.default(a)))
      temp[[1]] = NULL
    }
  }
  comps
}
