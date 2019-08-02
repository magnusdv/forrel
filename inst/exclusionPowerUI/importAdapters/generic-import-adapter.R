#' Get the index of a marker in a pedigree optionally creating it if it does not
#' exist.
#'
#' @param ped a `ped` object
#' @param markerName the name of a marker to find in the pedigree
#'
#' @return a 2-tuple in the form of a list containing the posibly updated ped
#'   and the index of the sought marker
#'
#' @example
#' p = nuclearPed(1)
#'
#' res = getOrAttachMarker(p, 'TH01')
#' p = res$ped
#' index = res$index
#' # marker TH01 is now guaranteed to be at ped$markerdata[[index]]
getOrAttachMarker = function(ped, markerName) {
  foundMarkers = tryCatch(
    whichMarkers(ped, markerName),
    error = function(e) NULL
  )

  if (length(foundMarkers) == 0) {
    m = marker(ped, name = markerName)
    ped = addMarkers(ped, m)

    list(ped = ped, index = whichMarkers(ped, c(markerName))[1])
  } else if (length(foundMarkers) == 1) {
    list(ped = ped, index = foundMarkers[1])
  } else {
    stop("Refusing to attach new markers to pedigree with duplicate markers.")
  }
}

#' Attach allele denomination and frequency data to a pedigree
#'
#' Data is read from a given dataframe, or from a file. If a dataframe is given,
#' the expected format is one marker per column and one allele per row, with
#' colnames on the first row and rownames on the first column. Variations of
#' this format can be achieved by modifying the parameters in ..., which are
#' passed directly to [utils::read.table()]. See examples for more details.
#'
#' If there are markers already in the pedigree, they are updated with new
#' frequency data if a matching marker is found in the provided dataset. If
#' there are duplicate markers in the given pedigree, this function will fail.
#'
#' @param ped a `ped` object
#' @param markers a list of markers to process
#' @param df a dataframe to read data from, defaults to NULL
#' @param Xchrom a list of marker names that are on the X chromosome
#' @param ... optional parameters to be passed to [utils::read.table()]. At
#'   least one of `file` or `text` is required if df is NULL.
#'
#' @return the `ped` object with attached genetic data
#'
#' @author Elias Hernandis <eliashernandis@gmail.com>
#'
#' @seealso [pedtools::readPed()]
#'
#' @examples
#' # read allele data from a CSV which we have stored in a string
#' csv = ",\"M1\",\"M2\"\n1,0.5,NA\n2,0.5,0.2\n3,NA,0.8"
#' ped = nuclearPed(1)
#' ped = attachAlleleFrequenciesToPedigree(ped, text = csv, sep = ',', header = TRUE, row.names = 1)
#'
#' # read data first to a dataframe (allows for more advanced pre-processing)
#' # and then attach it to the pedigree
#' df = read.table(text = csv, sep = ',', header = TRUE, row.names = 1)
#' ped = nuclearPed(1)
#' ped = attachAlleleFrequenciesToPedigree(ped, df = df)
#'
#' # read only some markers (leaving existing ones untouched)
#' ped = attachAlleleFrequenciesToPedigree(ped, markers = c('M1'), df = df)

attachAlleleFrequenciesToPedigree = function(ped, markers = NULL, df = NULL, Xchrom = NULL, ...) {
  if (is.pedList(ped)) {
    return(lapply(ped, function(x) {
      attachAlleleFrequenciesToPedigree(x, markers, df, Xchrom, ...)
    }))
  }

  if (is.null(df)) {
    df = read.table(...)
  }

  if (is.null(markers)) {
    markers = colnames(df)
  }

  for (markerName in markers) {
    res = getOrAttachMarker(ped, markerName)
    ped = res$ped
    index = res$index

    # attach allele denominations and frequency data to the marker
    als = rownames(df[!is.na(df[markerName]), , drop = FALSE])

    # get frequency data
    freqs = df[!is.na(df[markerName]), markerName]

    attr(ped$markerdata[[index]], 'alleles') = als
    attr(ped$markerdata[[index]], 'afreq') = freqs

    if (markerName %in% Xchrom) {
      attr(ped$markerdata[[index]], 'chrom') = 23
    }
  }

  ped
}

attachGenotypeToPedigree = function(ped, markers = NULL, df = NULL, ...) {
  if (is.pedList(ped)) {
    return(lapply(ped, function(x) {
      attachGenotypeToPedigree(x, markers, df, ...)
    }))
  }

  if (is.null(df)) {
    df = read.table(...)
  }

  if (is.null(markers)) {
    markers = unique(df[,2])
  }

  for (markerName in markers) {
    # find df rows concerning markerName
    relevant = df[df[,2] == markerName,]
    for (person in relevant[,1]) {
      if (person %in% labels(ped)) {
        alleles = c(relevant[relevant[,1] == person, 3],
                    relevant[relevant[,1] == person, 4])
        ped = setAlleles(ped,
                         ids = person,
                         markers = markerName,
                         alleles = alleles)
      }
    }
  }

  ped
}
