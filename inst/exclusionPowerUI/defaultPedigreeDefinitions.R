# List of pedigrees that are shown in the claim and true pedigree dropdowns. The
# list values correspond to those supported by the [pedigreeFromUI()] function.
# Therefore, to add other default pedigrees to the list, you need to add them in
# that function in addition to here.
defaultPedigrees = c("Nuclear family (1 son)" = "nucPed-1s",
                     "Nuclear family (1 daughter)" = "nucPed-1d",
                     "Unrelated" = "unrelated",
                     "Custom (from .ped file)" = "pedfile")


#' Create an unrelated true pedigree with the individuals available from
#' genotyping from the claim pedigree.
#'
#' @param claimPed the claim pedigree as a `ped` object
#' @param ids the list of individuals available for genotyping as a list of
#'   strings
#'
#' @return a list of ped objects, each containing an individual from the id list
#'
#' @note the `claimPed` param is required to get the sexes for the members of
#'   the unrelated pedigree
#'
#' @references See
#'   https://github.com/magnusdv/forrel/wiki/Exclusion-Power-GUI-how-to-define-an-unrelated-pedigree-as-the-true-pedigree
#'
#' @author Elias Hernandis <eliashernandis@gmail.com>
unrelatedPedFromClaimPed = function(claimPed, ids) {
  if (length(ids) == 0) return(list())

  sexes = getSex(claimPed, ids)
  tuples = cbind(ids, sexes)
  return(apply(tuples, 1, function(row) {pedtools::singleton(row[1], sex = row[2])}))
}

#' Create a `ped` object from the selected item in the pedigree dropdown.
#'
#' @param pedigreeID an ID defined in the values of the `defaultPedigrees` list
#'   defined avobe
#' @param pedfile a pedigree file, ignored if `pedigreeID` is not `pedfile`
#'
#' @note Unrelated pedigrees are handled separately by the
#'   [unrelatedPedFromClaimPed()] function.
#'
#' @return a `ped` object
#'
#' @author Elias Hernandis <eliashernandis@gmail.com>
pedigreeFromUI = function(pedigreeID, pedfile = NULL) {
  # Nuclear family with one son
  if (pedigreeID == 'nucPed-1s') {
    return(nuclearPed(1, father = "Father", mother = "Mother", children = c("Son")))
  }
  # Nuclear family with one daughter
  else if (pedigreeID == 'nucPed-1d') {
    return (nuclearPed(1, sex = 2, father = "Father", mother = "Mother", children = c("Daughter")))
  }
  # From a pedigree file
  else if (pedigreeID == 'pedfile') {
    if (is.null(pedfile)) {
      return()
    }
    return(as.ped(read.table(pedfile$datapath)))
  }
}
