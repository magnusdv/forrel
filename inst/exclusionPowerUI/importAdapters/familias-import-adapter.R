attachAlleleFrequenciesToPedigree.familias = function(ped, path, markers = NULL, Xchrom = NULL, ...) {
  if (is.pedList(ped)) {
    return(lapply(ped, function(x) {
      attachAlleleFrequenciesToPedigree.familias(x, path, markers, Xchrom, ...)
    }))
  }

  conn = file(description = path, open = 'r')
  lines = readLines(conn)

  i = 1
  while (i <= length(lines)) {
    markerName = lines[i]
    i = i + 1

    alleles = c()
    afreq = c()

    while (i <= length(lines)) {
      if (lines[i] == '') break

      unpacked = unlist(strsplit(lines[i], '\t', fixed = TRUE))
      al = unpacked[1]
      freq = as.double(unpacked[2])

      alleles = c(alleles, al)
      afreq = c(afreq, freq)
      i = i + 1
    }

    i = i + 1

    if (is.null(markers) || markerName %in% markers) {
      res = getOrAttachMarker(ped, markerName)
      ped = res$ped
      index = res$index

      attr(ped$markerdata[[index]], 'alleles') = alleles
      attr(ped$markerdata[[index]], 'afreq') = afreq

      if (markerName %in% Xchrom) {
        attr(ped$markerdata[[index]], 'chrom') = 23
      }
    }
  }

  close(conn)

  ped
}

attachGenotypeToPedigree.familias = function(ped, markers = NULL, df = NULL, ...) {
  if (is.pedList(ped)) {
    return(lapply(ped, function(x) {
      attachGenotypeToPedigree.familias(x, markers, df, ...)
    }))
  }

  if (is.null(df)) {
    df = read.table(...)
  }

  if (is.null(markers)) {
    headers = colnames(df)[2:ncol(df)]
    print(headers)
  }

  ped
}
