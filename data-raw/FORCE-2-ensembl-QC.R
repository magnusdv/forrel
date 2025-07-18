library(tidyverse)

# Load downloaded Ensembl data
ens = readRDS("data-raw/ensemblFORCE.rds")

# Load FORCE data - without Y and 1 missing SNP
xraw = readxl::read_xlsx("data-raw/FORCE-data.xlsx",
                         skip = 2,
                         n_max = 5422,  # skip bad ones at the bottom
                         col_types = "text") |>
  filter(SNP %in% names(ens)) |>
  print()

# Extract population data
POPDATA = map(ens, ~as_tibble(.$populations))


# QC ----------------------------------------------------------------------

getPop = function(pop) {
  map(POPDATA, \(x) filter(x, population == pop) |>
        select(allele, frequency) |>
        arrange(-frequency) |>
        deframe())
}

kg = getPop("1000GENOMES:phase_3:ALL")
gnom = getPop("gnomADg:ALL")
lengths(kg) |> table()
lengths(gnom) |> table()

# gnomad has many triallelic, but these are tiny
map_dbl(gnom, \(x) x[3]) |> summary()

# New version of gnom
gnom = map(gnom, \(x) x[x > 0.01])
lengths(gnom) |> table() # better

# Utility
pst = function(v, v2 = NULL) {
  if(is.null(v2)) paste(v, collapse = "/")
  else paste(v, v2, sep = ":", collapse = "/")
}

# Collect various info for QC controls
als = map_dfr(names(kg), \(x) {
  kx = kg[[x]]
  gx = gnom[[x]]
  tibble(snp = x,
         n_kg = length(kx),
         max_kg = if(n_kg > 0) max(kx) else NA,
         kg = pst(names(kx), round(kx,3)),
         equal = setequal(names(kx), names(gx)),
         n_gn = length(gx),
         max_gn = if(n_gn > 0) max(gx) else NA,
         gnom = pst(names(gx), round(gx,3))
  )})

# Inspect bad ones:
filter(als, n_kg < 2 | max_kg > 0.9 | !equal) |> arrange(n_kg) |> print(n = Inf)

# Conclusions:
# * 12 where 1000g has < 2 alleles: Use gnomADg for all of these
# * Many where 1000g has maxf > 0.9: All look legit and can be used
# * Remaining where allele sets differ: 1000g look ok.

