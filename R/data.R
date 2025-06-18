#' Norwegian STR frequencies
#'
#' A database of Norwegian allele frequencies for 35 STR markers.
#'
#' @format A list of length 35. Each entry is a numerical vector summing to 1,
#'   named with allele labels.
#'
#'   The following markers are included:
#'
#'   * `D3S1358`: 12 alleles
#'
#'   * `TH01`: 10 alleles
#'
#'   * `D21S11`: 26 alleles
#'
#'   * `D18S51 `: 23 alleles
#'
#'   * `PENTA_E`: 21 alleles
#'
#'   * `D5S818`: 9 alleles
#'
#'   * `D13S317`: 9 alleles
#'
#'   * `D7S820`: 19 alleles
#'
#'   * `D16S539`: 9 alleles
#'
#'   * `CSF1PO`: 11 alleles
#'
#'   * `PENTA_D`: 24 alleles
#'
#'   * `VWA`: 12 alleles
#'
#'   * `D8S1179`: 12 alleles
#'
#'   * `TPOX`: 9 alleles
#'
#'   * `FGA`: 25 alleles
#'
#'   * `D19S433`: 17 alleles
#'
#'   * `D2S1338`: 13 alleles
#'
#'   * `D10S1248`: 9 alleles
#'
#'   * `D1S1656`: 17 alleles
#'
#'   * `D22S1045`: 9 alleles
#'
#'   * `D2S441`: 13 alleles
#'
#'   * `D12S391`: 23 alleles
#'
#'   * `SE33`: 55 alleles
#'
#'   * `D7S1517`: 11 alleles
#'
#'   * `D3S1744`: 8 alleles
#'
#'   * `D2S1360`: 10 alleles
#'
#'   * `D6S474`: 6 alleles
#'
#'   * `D4S2366`: 7 alleles
#'
#'   * `D8S1132`: 12 alleles
#'
#'   * `D5S2500`: 8 alleles
#'
#'   * `D21S2055`: 18 alleles
#'
#'   * `D10S2325`: 10 alleles
#'
#'   * `D17S906`: 78 alleles
#'
#'   * `APOAI1`: 41 alleles
#'
#'   * `D11S554`: 51 alleles
#'
#' @source Dupuy et al. (2013): *Frequency data for 35 autosomal STR markers in
#'   a Norwegian, an East African, an East Asian and Middle Asian population and
#'   simulation of adequate database size*. Forensic Science International:
#'   Genetics Supplement Series, Volume 4 (1).
#'
"NorwegianFrequencies"


#' FORCE panel kinship SNPs
#'
#' A data frame describing (a subset of) the FORCE panel of SNPs designed for
#' applications in forensic genetics (Tillmar et al., 2021). The subset included
#' here are the SNPs recommended for kinship analysis. As the original
#' publication did not include allele frequencies, these were downloaded from
#' Ensembl via the biomaRt package. 15 markers were removed as frequency
#' information could not be retrieved.
#'
#' To attach the FORCE markers to a pedigree, use [pedtools::setSNPs()] (see
#' Examples).
#'
#'
#' @format A data frame with 3915 rows and 6 columns:
#'
#'   * `CHROM`: Chromosome (1-22)
#'
#'   * `MARKER`: Marker name (rs number)
#'
#'   * `MB`: Physical position in megabases (build GRCh38)
#'
#'   * `A1`: First allele
#'
#'   * `A2`: Second allele
#'
#'   * `FREQ1`: Allele frequency of `A1`
#'
#' @source Tillmar et al. The FORCE Panel: An All-in-One SNP Marker Set for
#'   Confirming Investigative Genetic Genealogy Leads and for General Forensic
#'   Applications. Genes. (2021)
#'
#' @examples
#' x = setSNPs(nuclearPed(), snpData = FORCE)
#' summary(x)
#'
#' getMap(x, markers = 1:5)
#' getFreqDatabase(x, markers = 1:5)
#'
"FORCE"


#' FORCE panel kinship SNPs (X-chromosomal)
#'
#' A data frame describing the 246 X-chromosomal SNPs included on the FORCE
#' panel for forensic genetics (Tillmar et al., 2021). As the original
#' publication did not include allele frequencies, these were downloaded from
#' Ensembl via the biomaRt package. For 9 markers the frequencies were obtained
#' manually from dbSNP.
#'
#' To attach the XFORCE markers to a pedigree, use [pedtools::setSNPs()] (see
#' Examples).
#'
#' @format A data frame with 246 rows and 6 columns:
#'
#'   * `CHROM`: Chromosome (all equal to "X")
#'   * `MARKER`: Marker name (rs number)
#'   * `MB`: Physical position in megabases (build GRCh37)
#'   * `A1`: First allele
#'   * `A2`: Second allele
#'   * `FREQ1`: Allele frequency of `A1`
#'
#' @source Tillmar et al. The FORCE Panel: An All-in-One SNP Marker Set for
#'   Confirming Investigative Genetic Genealogy Leads and for General Forensic
#'   Applications. Genes. (2021)
#'
#' @examples
#' x = setSNPs(nuclearPed(), snpData = XFORCE)
#' summary(x)
#'
#' getMap(x, markers = 1:5)
#' getFreqDatabase(x, markers = 1:5)
#'
"XFORCE"
