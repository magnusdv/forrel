library(forrel)

source('defaultPedigreeDefinitions.R')

source('advancedTableFileLoader.R')

source('tabularDataPreview.R')

source('importAdapters/generic-import-adapter.R')

getMarkerNames = function(ped) {
  if (is.pedList(ped)) {
    return(getMarkerNames(ped[[1]]))
  }

  as.vector(lapply(ped$markerdata, function(item) attr(item, 'name')), mode = 'character')
}
