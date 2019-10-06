library(forrel)

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

#' Get the names of the markers in a pedigree as a list of strings
#'
#' @param ped a `ped` object, or a list of such
#' @return a list of marker names as strings
getMarkerNames = function(ped) {
  if (is.pedList(ped)) {
    return(getMarkerNames(ped[[1]]))
  }

  as.vector(lapply(ped$markerdata, function(item) attr(item, 'name')), mode = 'character')
}

#' Returns the union over all markers in a pedigree of the possible alleles on
#' each marker.
#'
#' @param ped a `ped` object or a list of such
#' @return a list of strings
allAlleles = function(ped) {
  if (is.pedList(ped)) {
    return(allAlleles(ped[[1]]))
  }

  list = c()
  for (marker in ped$markerdata) {
    list = union(list, attr(marker, 'alleles'))
  }

  sort(list)
}

#' Extracts allele frequencies and denominations from a pedigree and returns
#' them in a matrix-like format, where there is one column per marker and one
#' row per allele. Alleles which do not occur on a marker have frequency 0 or
#' NA.
#'
#' @param ped a `ped` object, or a list of such. If this parameter is a list,
#'   the data is extracted from the first pedigree only.
#' @return a [data.frame] with one marker per column and one allele per row.
getTabularFrequencyDb = function(ped) {
  if (is.pedList(ped)) {
    return(getTabularFrequencyDb(ped[[1]]))
  }

  df = data.frame('Allele' = allAlleles(ped))

  for (marker in ped$markerdata) {
    colI = vapply(df[[1]], function(allele) {
      als = attr(marker, 'alleles')
      afreq = attr(marker, 'afreq')

      if (allele %in% als) {
        afreq[match(allele, als)]
      } else {
        NA
      }
    }, 1)

    df = cbind(df, colI)
  }

  colnames(df) = c('Alleles', getMarkerNames(ped))

  as.data.frame(df[2:ncol(df)], row.names = as.vector(df[[1]], mode = 'character'))
}

#' Returns true if any of the markers in pedigree is defined for the given ID.
#'
#' @param ped a `ped` object
#' @param id a string identifying a member of the pedigree
#' @return TRUE if the individual is genotyped, FALSE otherwise.
isGenotyped = function(ped, id) {
  any(!is.na(getAlleles(ped, ids = c(id))))
}

#' Returns a list of the individuals in the given pedigree that are genotyped as
#' defined by the function [isGenotyped()], i.e. the list of individuals for
#' which the function [isGenotyped()] returns TRUE.
#'
#' @param ped a `ped` object, or a list of such.
#' @return a list of ids (strings) of the genotyped individuals.
getGenotypedIds = function(ped) {
  if (is.pedList(ped)) {
    g = c()
    for (x in ped) {
      g = c(g, getGenotypedIds(x))
    }

    return(g)
  }

  labels(ped)[as.vector(lapply(labels(ped), function(id) { isGenotyped(ped, id) }), mode = 'logical')]
}


source('defaultPedigreeDefinitions.R')
source('uiModules/advancedTableFileLoader.R')
source('uiModules/tabularDataPreview.R')
source('uiModules/markerSettingsTable.R')

source('importAdapters/generic-import-adapter.R')
source('importAdapters/familias-import-adapter.R')
