source('defaultPedigreeDefinitions.R')

source('advancedTableFileLoader.R')

source('importAdapters/generic-import-adapter.R')

getMarkerNames = function(ped) {
  as.vector(lapply(ped$markerdata, function(item) attr(item, 'name')), mode = 'character')
}
