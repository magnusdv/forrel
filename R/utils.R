stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a single number, with optional range constraints
is_number = function(x, minimum = NA, maximum = NA) {
  isTRUE(length(x) == 1 &&
           is.numeric(x) &&
           (is.na(minimum) || x >= minimum) &&
           (is.na(maximum) || x <= maximum))
}

pluralise = function(noun, n) {
  if(n == 1) noun else sprintf("%ss", noun)
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}

#random 0/1 vector of length n.
.rand01 = function(n) {
  sample.int(2, size = n, replace = TRUE) - 1L
}


# Faster version of t.default(combn(n, 2, simplify = T))
.comb2 = function(n, vec = length(n) > 1){
  if(vec) {
    v = n
    n = length(v)
  }
  if (n < 2)
    return(matrix(nrow = 0, ncol = 2))

  x = rep.int(seq_len(n - 1), (n - 1):1)
  o = c(0, cumsum((n-2):1))
  y = seq_along(x) + 1 - o[x]

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

# Test if genotypes are consistent with ped
# (A better, but slower, alternative to `mendelianCheck()`)
consistentMarkers = function(x, markers = seq_len(nMarkers(x))) {

  # `marker` may be numeric, character or logical
  x = selectMarkers(x, markers)
  nMark = if(is.logical(markers)) sum(markers) else length(markers)

  # Compute likelihoods with no mutation model
  liks = vapply(seq_len(nMark), function(i) {
    mutmod(x, i) = NULL
    pedprobr::likelihood(x, i)
  }, FUN.VALUE = 0)

  # Return TRUE if likelihood is nonzero
  liks > 0
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
