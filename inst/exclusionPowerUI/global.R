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

getMarkerNames = function(ped) {
  if (is.pedList(ped)) {
    return(getMarkerNames(ped[[1]]))
  }

  as.vector(lapply(ped$markerdata, function(item) attr(item, 'name')), mode = 'character')
}

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



source('defaultPedigreeDefinitions.R')
source('advancedTableFileLoader.R')
source('tabularDataPreview.R')

source('importAdapters/generic-import-adapter.R')
source('importAdapters/familias-import-adapter.R')
