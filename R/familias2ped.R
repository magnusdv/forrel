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
#' @references Windows Familias is freely available from <http://familias.name>.
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

  p = data.frame(id = id, fid = 0, mid = 0, sex = 1, stringsAsFactors = FALSE)
  p$fid[findex > 0] = id[findex]
  p$mid[mindex > 0] = id[mindex]
  p$sex[familiasped$sex == "female"] = 2

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
    addedParents = matrix("0", nrow = nFath + nMoth, ncol = ncol(datamatrix))
    allelematrix = rbind(allelematrix, addedParents)

    p = cbind(p, allelematrix, stringsAsFactors = FALSE)
  }

  ### Part 3: marker annotations

  annotations = readFamiliasLoci(loci)

  ### Create ped object

  pedtools::as.ped(p, locus_annotations = annotations)
}



#' @rdname Familias2ped
#' @importFrom pedmut mutationMatrix mutationModel validateMutationModel
#' @export
readFamiliasLoci = function(loci) {
  if (is.null(loci))
    return(NULL)
  if (class(loci) == "FamiliasLocus")
    loci = list(loci)

  lapply(loci, function(a) {
    als = names(a$alleles)
    afreq = as.numeric(a$alleles)

    malemut = a$maleMutationMatrix
    femalemut = a$femaleMutationMatrix

    if (all(diag(malemut) == 1))
      malemut = NULL
    else
      malemut = mutationMatrix("custom", matrix = malemut)

    if (all(diag(femalemut) == 1))
      femalemut = NULL
    else
      femalemut = mutationMatrix("custom", matrix = femalemut)

    if (is.null(malemut) && is.null(femalemut))
      mutmod = NULL
    else {
      mutmod = mutationModel(list(female = femalemut, male = malemut))
      validateMutationModel(mutmod, alleles = als)
    }

    list(name = a$locusname, alleles = als, afreq = afreq, mutmod = mutmod)
  })
}
