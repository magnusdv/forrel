library(tidyverse)

# Load downloaded Ensembl data
ens = readRDS("data-raw/ensemblFORCE.rds")

# Load FORCE data - without Y and 1 missing SNP
x = readxl::read_xlsx("data-raw/FORCE-data.xlsx",
                         skip = 2,
                         n_max = 5422,  # skip bad ones at the bottom
                         col_types = "text") |>
  filter(SNP %in% names(ens)) |>
  print()

# Extract population data
POPDATA = map(ens, ~as_tibble(.$populations))



# Create datasets ---------------------------------------------------------


# Based on QC conclusions we use the following freqs:
# * 1000G:ALL unless missing, otherwise gnomADg:ALL

FORCE_FREQ = map_dfr(POPDATA, .id = "MARKER", \(rs) {
  dat = filter(rs, population == "1000GENOMES:phase_3:ALL")
  if(nrow(dat) != 2)
    dat = filter(rs, population == "gnomADg:ALL")

  if(nrow(dat) < 2) stop("Missing data!") # should not happen

  dat = arrange(dat, -frequency)

  tibble(A1 = dat$allele[1], A2 = dat$allele[2], FREQ1 = dat$frequency[1],
         SOURCE = dat$population[1])
}) |> print()


# Merge with position data
FORCE_ALL = x |>
  dplyr::filter(Type %in% c("Kinship", "Kinship_iiSNP", "XSNP")) |>
  dplyr::rename(
    CHROM = Chromosome,
    MARKER = SNP,
    POS = "GRCh38 Region",
  ) |>
  left_join(FORCE_FREQ, by = "MARKER") |>
  dplyr::mutate(
    MB = as.integer(POS)/1e6,
    FREQ1 = round(FREQ1, 3)
  ) |>
  dplyr::select(CHROM, MARKER, MB, A1, A2, FREQ1) |>
  print()


table(fct_inorder(FORCE_ALL$CHROM)) # 1-22, X

# Separate datasets for autosomal and X
FORCE  = FORCE_ALL |> filter(CHROM %in% 1:22) |> print()
XFORCE = FORCE_ALL |> filter(CHROM == "X")    |> print()

# Save datasets
usethis::use_data(FORCE,  overwrite = TRUE)
usethis::use_data(XFORCE, overwrite = TRUE)


# POP-SPECIFIC FREQS ------------------------------------------------------


# Function for adding column with pop-specific frequencies
setPopFreqs = function(x, popdata = POPDATA) {
  a1 = x$A1
  global = x$FREQ1
  names(a1) = names(global) = x$MARKER

  popfr = vapply(x$MARKER, function(r) {
    a = a1[r]
    b = popdata[[r]]
    bb = b[b$population == pop, , drop = FALSE]
    if(!nrow(bb))
      return(NA_real_)

    bb$frequency[match(a1[r], bb$allele)]
  }, numeric(1))

  x[[sub(".*:", "", pop)]] = popfr
  x
}

