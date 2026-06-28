stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a single number, with optional range constraints
isNumber = function(x, minimum = NA, maximum = NA) {
  isTRUE(length(x) == 1 &&
           is.numeric(x) &&
           (is.na(minimum) || x >= minimum) &&
           (is.na(maximum) || x <= maximum))
}

# Faster alternative to suppressWarnings(as.numeric(v))
# NB: doesn't catch scientific notation 1e-4.
asNum = function(v) {
  if(is.numeric(v)) return(v)
  num = grep("[^-0-9.]", v, invert = TRUE)
  u = rep(NA_real_, length(v))
  u[num] = as.numeric(v[num])
  u
}

# Faster alternative to suppressWarnings(as.integer(v))
# NB: doesn't catch scientific notation 1e+12.
asInt = function(v) {
  if(is.numeric(v)) return(as.integer(v))
  num = grep("[^-0-9]", v, invert = TRUE)
  u = rep(NA_integer_, length(v))
  u[num] = as.integer(v[num])
  u
}

# round + toString
rst = function(v, digits = 3)
  toString(round(v, digits))

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

`%notin%` = function(x, table)
  match(x, table, nomatch = 0L) == 0L

pluralise = function(noun = "", n) {
  if(missing(n)) return(pluralise("", noun))
  if(n == 1) noun else sprintf("%ss", noun)
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}

.mysetequal = function(x, y) {
  !anyNA(match(x, y)) && !anyNA(match(y, x))
}

#random 0/1 vector of length n.
.rand01 = function(n) {
  sample.int(2, size = n, replace = TRUE) - 1L
}

.setnames = function(x, nms) {
  names(x) = nms
  x
}

.miraiWorkers = function() {
  if(!mirai::daemons_set())
    return(0L)
  mirai::info()[["connections"]] %||% 0L
}

ftime = function(st, digits = 3)
  format(Sys.time() - st, digits = digits)


# Faster version of t.default(combn(n, 2, simplify = T))
.comb2 = function(n, vec = length(n) > 1){
  if(vec) {
    v = n
    n = length(v)
  }
  if(n < 2)
    return(matrix(nrow = 0, ncol = 2))

  x = rep.int(seq_len(n - 1), (n - 1):1)
  y = sequence.default((n - 1):1, 2:n)

  if(vec)
    cbind(v[x], v[y], deparse.level = 0)
  else
    cbind(x, y, deparse.level = 0)
}

isEP = function(x) {
  inherits(x, "EPresult") || inherits(x, "mpEP")
}

isIP = function(x) {
  inherits(x, "LRpowerResult") || inherits(x, "mpIP")
}


# Disable mutations
disableMutationModels = function(x, disable, verbose = FALSE) {

  if(isFALSE(disable) || is.null(disable))
    return(x)

  # Which of the markers allow mutations?
  hasMut = allowsMutations(x)

  # Return early if no markers has mutation models
  if(!any(hasMut))
    return(x)

  if(isTRUE(disable))
    disable = which(hasMut)
  else if(identical(disable, NA)) # Disable for consistent markers
    disable = which(hasMut & consistentMarkers(x, hasMut))
  else # if numeric or character
    disable = whichMarkers(x, disable)

  # Disable
  if(length(disable)) {
    if(verbose)
      message("Disabling mutations for markers: ", toString(disable))
    mutmod(x, disable) = NULL
  }

  # Return the modified object
  x
}


# TODO: Move to pedtools?
fixAllelesAndFreqs = function(alleles = NULL, afreq = NULL,
                              observed = NULL, NAstrings = c(0, "", NA, "-")) {

  if(!is.null(alleles) && !is.null(names(afreq)))
    stop2("Argument `alleles` should not be used when `afreq` has names")
  if(is.null(alleles) && !is.null(afreq) && is.null(names(afreq)))
    stop2("When `alleles` is NULL, `afreq` must be named")

  # If alleles are NULL, take from afreq names, otherwise from supplied genos
  als = alleles %||% names(afreq) %||% .mysetdiff(observed, NAstrings)
  if(length(als) == 0)
    als = 1:2

  ### Frequencies
  afreq = afreq %||% {rep_len(1, length(als))/length(als)}
  names(afreq) = names(afreq) %||% als

  # Sort alleles and frequencies (numerical sorting if appropriate)
  numAls = asNum(als)
  ord = if(!anyNA(numAls)) ord = order(numAls) else order(als)

  # Return ordered, named frequencies
  afreq[ord]
}


# TODO: Remove the following

#' Add points to the IBD triangle
#'
#' This function is re-exported from the `ribd` package. For documentation see
#' [ribd::showInTriangle()].
#'
#' @importFrom ribd showInTriangle
#' @name showInTriangle
#' @export
NULL
