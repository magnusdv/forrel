#' Import/export from Familias
#'
#' @description
#' Functions for reading .fam files associated with the Familias
#' software for forensic kinship computations.
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

