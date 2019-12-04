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
  sample.int(2, size = n, replace = T) - 1L
}

# Equivalent to t.default(combn(n, 2)), but ~6 times faster.
.comb2 = function(n) {
  if (n < 2)
    return(matrix(nrow = 0, ncol = 2))
  v1 = rep.int(seq_len(n - 1), (n - 1):1)
  v2 = NULL
  for (i in 2:n) v2 = c(v2, i:n)
  cbind(v1, v2, deparse.level = 0)
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
