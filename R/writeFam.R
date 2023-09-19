#' Convert `ped` objects to `Familias` format
#'
#' This function produces a .fam file readable by the `Familias` software
#' (Egeland, Mostad et al., 2000), containing all input pedigrees and their
#' marker data. The option `openFam = TRUE` calls `openFamilias()` to open a
#' fresh `Familias` session with the produced file pre-loaded.
#'
#' The following parameters are applied by default, but may be adjusted with the
#' `params` argument:
#'
#' * `dbName = "unknown"`
#' * `dbSize = 1000`
#' * `dropout = 0`
#' * `maf = 0`
#' * `theta = 0`
#'
#' The `params` argument should be a list similar to the `params` slot produced
#' by [readFam()] with `includeParams = TRUE`. Single entries are recycled if
#' needed. If `params` contains a vector `dropout` with dropout probabilities
#' for certain pedigree members, it is converted into corresponding
#' `dropoutConsider` and `dropoutValue` vectors (see Examples).
#'
#' @param ... One or several pedigrees. Each argument should be either a single
#'   `ped` object or a list of such. If the pedigrees are unnamed, they are
#'   assigned names "Ped 1", "Ped 2", etc.
#' @param famfile The name or path to the output file to be written. The
#'   extension ".fam" is added if missing.
#' @param params A list of further parameters; see [readFam()] for valid
#'   entries. See also Details for default values.
#' @param dbOnly A logical. If TRUE, no pedigree information is included; only
#'   the frequency database.
#' @param openFam A logical. If TRUE, an attempt is made to open the produced
#'   `fam` file in an external `Familias` session. Only available on Windows
#'   systems with a working `Familias` installation.
#' @param FamiliasPath The path to the Familias executable. If empty, the
#'   following are tried in order: "Familias3.exe", "C:/Program Files
#'   (x86)/Familias3/Familias3.exe".
#' @param verbose A logical, by default TRUE.
#'
#' @return The filename is returned invisibly.
#'
#' @seealso [readFam()]
#'
#' @references Egeland, T., P. F. Mostad, et al. (2000). _Beyond traditional
#'   paternity and identification cases. Selecting the most probable pedigree._
#'   Forensic Sci Int 110(1): 47-59.
#'
#' @examples
#'
#' # Create pedigree with 2 markers
#' x = nuclearPed() |>
#'   profileSim(markers = NorwegianFrequencies[1:2], seed = 1)
#'
#' # Write to .fam
#' tmp = writeFam(x, famfile = tempfile())
#'
#' # Read back in
#' y = readFam(tmp)
#'
#' stopifnot(identical(x, y))
#'
#'
#' ### With stepwise mutation model
#' x2 = setMutmod(x, model = "stepwise",
#'                rate = list(male = 0.001, female = 0.002),
#'                range = 0.1, rate2 = 0.0001)
#'
#' # Write and read
#' y2 = x2 |>
#'   writeFam(famfile = tempfile()) |>
#'   readFam()
#'
#' stopifnot(identical(x2, y2))
#'
#'
#' ### Read/write including detailed parameters
#' params = list(theta = 0.1, dbName = "myDB", dropout = c("3" = 0.01))
#' fam = writeFam(x2, famfile = tempfile(), params = params)
#'
#' dat = readFam(fam, includeParams = TRUE)
#'
#' # Pedigree is now in the `main` slot
#' stopifnot(identical(x2, dat$main))
#'
#' # The `dropout` parameter is converted to (and is equivalent to):
#' dat$params$dropoutConsider
#' dat$params$dropoutValue
#'
#'
#' ### Read/write frequency database
#'
#' # Write database as fam file
#' dbfam = writeFam(x2, famfile = tempfile(), dbOnly = TRUE)
#'
#' # Read back in: The result is a list of marker attributes
#' a = readFam(dbfam)
#'
#' # Attach to a pedigree and write to a new file
#' z = singleton(1) |> setMarkers(locusAttributes = a)
#' dbfam2 = writeFam(z, famfile = tempfile(), dbOnly = TRUE)
#'
#' stopifnot(identical(readLines(dbfam), readLines(dbfam2)))
#'
#' @export
writeFam = function(..., famfile = "ped.fam", params = NULL, dbOnly = FALSE,
                    openFam = FALSE, FamiliasPath = NULL, verbose = TRUE) {
  peds = list(...)
  if (length(peds) == 1)
    peds = peds[[1]]
  if (is.ped(peds))
    peds = list(peds)

  # Ensure each entry is a pedlist
  peds = lapply(peds, function(p) if(is.ped(p)) list(p) else p)
  if(!all(isped <- sapply(peds, is.pedList))) {
    idx = which(!isped)[1]
    stop2(sprintf("Argument %d is not a pedigree, but: '%s'", idx,
                  class(peds[[idx]])))
  }

  # All unique individual names
  LABS = unique.default(unlist(lapply(peds, labels)))
  nind = length(LABS)

  # All unique marker names
  MARKERS = unique.default(unlist(lapply(peds, name)))
  nmar = length(MARKERS)

  # Extra param: Dropout
  dropoutConsider = params$dropoutConsider %||% setNames(rep_len(FALSE, nind), LABS)
  dropoutValue    = params$dropoutValue %||% setNames(rep_len(0, nmar), MARKERS)

  # If shortcut "dropout" is used, overrule the others
  if(!is.null(dropoutInd <- params[["dropout"]])) {

    # Familias doesn't support individual dropout values
    if(length(udr <- unique.default(dropoutInd[dropoutInd > 0])) > 1)
      stop2("All nonzero dropout values must be equal: ", sort(udr))

    # Override dropoutValue
    dropoutValue[] = max(udr)

    # Override params$dropoutConsider
    if(!is.null(dnms <- names(dropoutInd))) {
      if(!all(dnms %in% LABS))
        stop2("Unknown ID in `dropout`: ", setdiff(dnms, LABS))
      dropoutConsider[] = FALSE
      dropoutConsider[dnms] = dropoutInd > 0
    }
    else if(length(dropoutInd) == 1 && is.numeric(dropoutInd))
      dropoutConsider[] = dropoutInd > 0
  }

  # Make sure only typed individuals have dropout (otherwise Familias crashes!)
  untyped = unlist(lapply(peds, untypedMembers))
  dropoutConsider[untyped] = FALSE

  # Extra parameter: Database size for each marker
  dbSize = params$dbSize %||% 1000
  if(length(dbSize) == 1)
    dbSize = setNames(rep_len(dbSize, nmar), MARKERS)

  # Extra parameter: Minor allele frequency for each marker
  maf = params$maf %||% 0
  if(length(maf) == 1)
    maf = setNames(rep_len(maf, nmar), MARKERS)

  # Open output connection
  if(!endsWith(famfile, ".fam"))
    famfile = paste0(famfile, ".fam")
  fam = file(famfile, "w", encoding = "UTF-8")

  isClosed = FALSE
  on.exit(if(!isClosed) close(fam), add = TRUE)

  # Quick utilities
  addline = function(...) cat(..., file = fam, sep = "\n")
  quo = function(s) sprintf('"%s"', s)

  # Preamble
  addline(quo("Output from Familias, version 3.2.8"),
          quo(sprintf("(Actually produced by R/pedsuite, %s)", format(Sys.Date(), "%d %b %Y"))),
          "3.2.8",
          '""',
          if(dbOnly) 0 else nind)

  # Individuals and genotypes ---------------------------------------------

  # Loop over indivs
  taken = rep(FALSE, length(LABS))
  names(taken) = LABS

  # Hack to skip everyone if `dbOnly = TRUE`
  if(dbOnly)
    taken[] = TRUE

  for(ped in peds) for(x in ped) for(id in labels(x)) {
    if(taken[id])
      next

    addline(quo(id),
            "#FALSE#",
            paste(-1, if(dropoutConsider[id] > 0) "(Consider dropouts)"),
            "#FALSE#",
            ifelse(getSex(x, id) == 1, "#TRUE#", "#FALSE#"))

    idx = internalID(x, id)

    # 1-based list of triples (a1, a2, midx)
    alslist = lapply(x$MARKERS, function(m) {
      a = m[idx, ]
      if(any(a > 0)) {
        mIdx = match(attr(m, "name"), MARKERS)
        c(a, mIdx) # else NULL
      }
    })
    alsvec = unlist(alslist)

    # number of markers for this indiv
    k = length(alsvec) / 3
    addline(k)
    if(k > 0)
      addline(alsvec - 1)  # convert to 0-based!

    taken[id] = TRUE
  }

  # Pedigrees ---------------------------------------------------------------

  # Skip all pedigrees if `dbOnly = TRUE`
  if(dbOnly) {
    npeds = 0
  } else {
    npeds = length(peds)
    pednames = names(peds) %||% paste("Ped", 1:npeds)
  }

  addline('"Known relations"', 0,0,0, npeds)

  # Loop over pedigrees
  for(i in seq_len(npeds)) {
    ped = peds[[i]]
    nrels = sum(lengths(lapply(ped, nonfounders))) * 2

    addline(i - 1, quo(pednames[i]), 0, 0, nrels)

    for(x in ped) {
      nonfou = nonfounders(x, internal = TRUE)
      if(!length(nonfou))
        next

      # Parent-child relationships: 1-based index within component
      relx = rbind(x$FIDX[nonfou], nonfou,
                   x$MIDX[nonfou], nonfou, deparse.level = 0)

      # Convert to 1-based index within complete LABS
      relxx = match(x$ID[as.vector(relx)], LABS)

      # Convert to 0-based
      addline(relxx - 1)
    }
  }

  # Database -------------------------------------------------------

  addline(sprintf("#FALSE# (#Databases: 1 ;Theta/Kinship/Fst: %g )", params$theta %||% 0),
          nmar,
          "#TRUE#",
          quo(params$dbName))

  takenMark = rep(FALSE, length(MARKERS))
  names(takenMark) = MARKERS

  for(ped in peds) for(x in ped) for(m in x$MARKERS) {
    attrs = attributes(m)
    mname = attrs$name
    if(takenMark[mname])
      next

    # Extract mutation parameters
    hasMut = !is.null(attrs$mutmod)

    if(hasMut)
      mut = pedmut::getParams(attrs$mutmod, format = 1,
                              params = c("model", "rate", "range", "rate2"))
    else
      mut = data.frame(model = c("equal", "equal"), rate = 0, range = 0, rate2 = 0)

    mod = match(mut$model, c("equal", "proportional", "stepwise"), nomatch = 1L)
    mod[mod == 3 & !is.na(mut$rate2)] = 5

    nals = length(attrs$alleles)

    # Write marker name and mutation parameters
    addline(quo(mname),
            mut$rate %na% 0,
            mod - 1,
            nals,
            mut$range %na% 0,
            mut$rate2 %na% 0)

    addline("#FALSE#", 0,
            sprintf("%d\t(DatabaseSize = %d , Dropout probability = %g , Minor allele frequency = %g )",
                    nals, dbSize[[mname]], dropoutValue[[mname]], maf[[mname]]))

    frvec = character(2*nals)
    frvec[2*(1:nals) - 1] = quo(attrs$alleles)
    frvec[2*(1:nals)] = attrs$afreq

    addline(frvec)
    takenMark[mname] = TRUE
  }

  close(fam)
  isClosed = TRUE
  pth = normalizePath(famfile)

  if(verbose)
    cat("Written to file:", pth, "\n")

  if(openFam) {
    if(verbose)
      cat("Trying to open in Familias\n")
    openFamilias(pth, FamiliasPath, verbose = verbose)
  }

  invisible(pth)
}

#' @rdname writeFam
#' @export
openFamilias = function(famfile = NULL, FamiliasPath = NULL, verbose = TRUE) {
  FamiliasPath = .checkFamilias(FamiliasPath)
  if(verbose)
    cat("Familias executable:", FamiliasPath, "\n")
  if(is.null(famfile))
    cmd = sprintf('"%s"', shQuote(FamiliasPath))
  else
    cmd = sprintf('"%s %s"', shQuote(FamiliasPath), shQuote(famfile))

  shell(cmd, wait = FALSE)
}

.checkFamilias = function(FamiliasPath = NULL) {
  if(Sys.info()["sysname"] != "Windows")
    stop2("Familias is only available on Windows systems.")
  if(is.null(FamiliasPath))
    FamiliasPath = c("Familias3.exe", "C:/Program Files (x86)/Familias3/Familias3.exe")
  for(p in FamiliasPath) {
    if(file.exists(p) && file.info(p)$exe != "no")
      return(p)
  }
  stop2("Could not find Familias executable: ", FamiliasPath)
}

`%na%` = function(x, y) {
  if(length(x) == 1) {
    if(is.na(x)) y else x
  }
  else {
    x[is.na(x)] = y
    x
  }
}

.internalAlleleMatrix = function(x) {
  mat = matrix(unlist(x$MARKERS),
               nrow = length(x$ID),
               dimnames = list(x$ID, NULL))
  mat
}

