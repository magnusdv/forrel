stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call.=FALSE))
  do.call(stop, a)
}

pluralise = function(noun, n) {
  if(n==1) noun else sprintf("%ss", noun)
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
  sample.int(2, size = n, replace = T) - 1
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
