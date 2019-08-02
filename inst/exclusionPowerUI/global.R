library(forrel)

source('defaultPedigreeDefinitions.R')

source('advancedTableFileLoader.R')

source('tabularDataPreview.R')

source('importAdapters/generic-import-adapter.R')
source('importAdapters/familias-import-adapter.R')

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
