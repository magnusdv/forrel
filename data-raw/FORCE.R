library(tidyverse)
library(biomaRt)

xraw = readxl::read_xlsx("data-raw/FORCE-data.xlsx",
                         skip = 2,
                         n_max = 5422,  # skip bad ones at the bottom
                         col_types = "text")

x = xraw %>%
  dplyr::filter(startsWith(Type, "Kinship")) %>%
  dplyr::select(CHROM = Chromosome, MARKER = SNP, POS = "GRCh38 Region",
         REF = "Reference Allele", ALT = "Alternate Allele") %>%
  dplyr::mutate(POS = as.integer(POS)) %>%
  print


mart = useEnsembl(biomart = "snps", dataset = "hsapiens_snp")
info = getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start", "minor_allele", "minor_allele_freq"),
  filters = c("snp_filter", "chr_name"),
  values = list(x$MARKER, 1:22),
  mart = mart
) %>% as_tibble() %>% print()

filter(info, minor_allele == "")
anyDuplicated(info$refsnp_id)

y = left_join(x, info, by = c(MARKER = "refsnp_id")) %>%
  dplyr::filter(minor_allele != "") %>%
  dplyr::mutate(
    MB = POS / 1e6,
    A1 = if_else(minor_allele == REF, ALT, REF),
    A2 = if_else(minor_allele == REF, REF, ALT),
    FREQ1 = 1 - minor_allele_freq
  ) # %>% print

# Checks
filter(y, POS != chrom_start)
filter(y, is.na(FREQ1))
filter(y, nchar(A1) !=1 | nchar(A2) !=1)

# Final table
FORCE = y %>%
  dplyr::select(CHROM, MARKER, MB, A1, A2, FREQ1) %>%
  as.data.frame()

print(head(FORCE))

usethis::use_data(FORCE, overwrite = TRUE)
