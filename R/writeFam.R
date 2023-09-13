#' Convert `ped` objects to `Familias` format
#'
#' This function produces a .fam file readable by the `Familias` software
#' (Egeland, Mostad et al., 2000), containing all input pedigrees and their
#' marker data. The option `openFam = TRUE` calls `openFamilias()` to open a
#' fresh `Familias` session with the produced file pre-loaded.
#'
#' @param ... One or several pedigrees. Each argument should be either a single
#'   `ped` object or a list of such. If the pedigrees are unnamed, they are
#'   assigned names "Ped 1", "Ped 2", etc.
#' @param famfile The name or path to the output file to be written. The
#'   extension ".fam" is added if missing.
#' @param theta A number between 0 and 1 inclusive, indicating a theta
#'   correction for the marker database. By default 0.
#' @param dropout A number between 0 and 1 inclusive, or a named vector of such
#'   numbers, indicating dropout probability. By default 0.
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
#' @references Egeland, T., P. F. Mostad, et al. (2000). _Beyond traditional
#'   paternity and identification cases. Selecting the most probable pedigree._
#'   Forensic Sci Int 110(1): 47-59.
#'
#' @examples
#' library(pedprobr)
#'
#' x = nuclearPed(father = "AF", mother = "MO", children = "CH") |>
#'   profileSim(N = 1, ids = c("AF", "CH"), seed = 111,
#'              markers = NorwegianFrequencies[1:2])
#'
#' # Write to .fam
#' tmp = writeFam(x, famfile = tempfile())
#'
#' # Read back in
#' y = readFam(tmp)
#'
#' # Verify that likelihoods agree
#' stopifnot(all.equal(likelihood(x),
#'                     likelihood(y)))
#'
#'
#' ### With stepwise mutation model
#' x2 = setMutmod(x, model = "stepwise",
#'                rate = list(male = 0.001, female = 0.002),
#'                range = 0.1, rate2 = 0.0001)
#'
#' y2 = x2 |> writeFam(famfile = tempfile()) |> readFam()
#'
#' stopifnot(all.equal(likelihood(x2), likelihood(y2)))
#'
#' @export
writeFam = function(..., famfile = "ped.fam", theta = 0, dropout = 0,
                    openFam = FALSE, FamiliasPath = NULL, verbose = TRUE) {
  peds = list(...)
  if (length(peds) == 1)
    peds = peds[[1]]
  if (is.ped(peds))
    peds = list(peds)

  npeds = length(peds)
  pednames = names(peds) %||% paste("Ped", 1:npeds)

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

  # Dropout
  if(!is.null(dnms <- names(dropout))) {
    if(!all(dnms %in% LABS))
      stop2("Unknown ID in `dropout`: ", setdiff(dnms, LABS))
    rest = rep_len(0, nind - length(dropout))
    names(rest) = setdiff(LABS, dnms)
    dropout = c(dropout, rest)[LABS] # sort to look nice
  }
  else {
    dropout = rep_len(dropout, nind)
    names(dropout) = LABS
  }

  dropoutValue = max(dropout)

  if(length(udr <- unique.default(dropout[dropout > 0])) > 1)
    stop2("All nonzero dropout values must be equal: ", sort(udr))

  # Make sure only typed individuals have dropout (otherwise Familias crashes!)
  untyped = unlist(lapply(peds, untypedMembers))
  dropout[untyped] = 0

  # Unique marker names
  MARKERS = unique.default(unlist(lapply(peds, name)))

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
          nind)

  # Individuals and genotypes ---------------------------------------------

  # Loop over indivs
  taken = rep(FALSE, length(LABS))
  names(taken) = LABS

  for(ped in peds) for(x in ped) for(id in labels(x)) {
    if(taken[id])
      next

    addline(quo(id),
            "#FALSE#",
            if(dropout[id] > 0) "-1 (Consider dropouts)" else "-1",
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

  addline('"Known relations"', 0,0,0, npeds)

  # Loop over pedigrees
  for(i in 1:npeds) {
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

  addline(sprintf("#FALSE# (#Databases: 1 ;Theta/Kinship/Fst: %g )", theta),
          length(MARKERS),
          "#TRUE#",
          quo("Unknown_db"))

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
            sprintf("%d\t(DatabaseSize = 1000 , Dropout probability = %f , Minor allele frequency = 0 )", nals, dropoutValue))

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

