#' Import/export from Familias
#'
#' @description
#' Functions for reading and writing .fam files associated with the Familias
#' software for forensic kinship computations, and for converting objects from
#' the R package `Familias` into the `pedsuite`.
#'
#' ***Deprecated*** These functions have been moved to a separate package,
#' `pedFamilias`, and will be removed from `forrel` in a future version.
#'
#' @param ... Arguments passed on to the respective `pedFamilias` function.
#'
#' @name familias
NULL


#' @rdname familias
#' @export
readFam = function(...) {
  cat("Deprecated, use `pedFamilias::readFam()` instead.\n")
  pedFamilias::readFam(...)
}


#' @rdname familias
#' @export
writeFam = function(...) {
  cat("Deprecated, use `pedFamilias::writeFam()` instead.\n")
  pedFamilias::writeFam(...)
}


#' @rdname familias
#' @export
openFamilias = function(...) {
  cat("Deprecated, use `pedFamilias::openFamilias()` instead.\n")
  pedFamilias::openFamilias(...)
}


#' @rdname familias
#' @export
Familias2ped = function(...) {
  cat("Deprecated, use `pedFamilias::Familias2ped()` instead.\n")
  pedFamilias::Familias2ped(...)
}



#' @rdname familias
#' @export
readFamiliasLoci = function(...) {
  cat("Deprecated, use `pedFamilias::readFamiliasLoci()` instead.\n")
  pedFamilias::readFamiliasLoci(...)
}
