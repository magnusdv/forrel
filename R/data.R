#' Norwegian STR frequencies
#'
#' A database of Norwegian allele frequencies for 35 STR markers. NOTE: An
#' updated database of these and other markers, are available in the `norSTR`
#' package. The `NorwegianFrequencies` object is kept here for backward
#' compatibility.
#'
#'
#' @format A list of length 35. Each entry is a numerical vector summing to 1,
#'   named with allele labels.
#'
#'   The following markers are included:
#'
#'   * `D3S1358`: 12 alleles
#'   * `TH01`: 10 alleles
#'   * `D21S11`: 26 alleles
#'   * `D18S51 `: 23 alleles
#'   * `PENTA_E`: 21 alleles
#'   * `D5S818`: 9 alleles
#'   * `D13S317`: 9 alleles
#'   * `D7S820`: 19 alleles
#'   * `D16S539`: 9 alleles
#'   * `CSF1PO`: 11 alleles
#'   * `PENTA_D`: 24 alleles
#'   * `VWA`: 12 alleles
#'   * `D8S1179`: 12 alleles
#'   * `TPOX`: 9 alleles
#'   * `FGA`: 25 alleles
#'   * `D19S433`: 17 alleles
#'   * `D2S1338`: 13 alleles
#'   * `D10S1248`: 9 alleles
#'   * `D1S1656`: 17 alleles
#'   * `D22S1045`: 9 alleles
#'   * `D2S441`: 13 alleles
#'   * `D12S391`: 23 alleles
#'   * `SE33`: 55 alleles
#'   * `D7S1517`: 11 alleles
#'   * `D3S1744`: 8 alleles
#'   * `D2S1360`: 10 alleles
#'   * `D6S474`: 6 alleles
#'   * `D4S2366`: 7 alleles
#'   * `D8S1132`: 12 alleles
#'   * `D5S2500`: 8 alleles
#'   * `D21S2055`: 18 alleles
#'   * `D10S2325`: 10 alleles
#'   * `D17S906`: 78 alleles
#'   * `APOAI1`: 41 alleles
#'   * `D11S554`: 51 alleles
#'
#' @source Dupuy et al. (2013): *Frequency data for 35 autosomal STR markers in
#'   a Norwegian, an East African, an East Asian and Middle Asian population and
#'   simulation of adequate database size*. Forensic Science International:
#'   Genetics Supplement Series, Volume 4 (1).
#'
"NorwegianFrequencies"


#' FORCE panel SNP data
#'
#' Data frames describing the FORCE panel of SNPs for forensic genetics (Tillmar
#' et al., 2021). We provide here two subsets of the complete panel: the
#' autosomal kinship SNPs: The autosomal kinship SNPs (`FORCE`, n = 3930) and
#' the X-chromosomal SNPs (`XFORCE`, n = 246). To attach the markers to a
#' pedigree, use [pedtools::setSNPs()] (see Examples).
#'
#' Allele frequencies were retrieved from Ensembl using the REST API, with the
#' population `1000GENOMES:phase_3:ALL` as primary source. For 9 SNPs where this
#' was unavailable, `gnomADg:ALL` was used instead. The SNP rs2323964 was
#' excluded due to lack of Ensembl support.
#'
#' The autosomal dataset (`FORCE`) was updated in version 1.8.1, adding 15
#' markers that were previously missing and revising some frequencies. The
#' previous version is available via `system.file("FORCE_old", package =
#' "forrel")`.
#'
#' For details, the code used to download and process the data is available in
#' the `data-raw` folder on GitHub:
#' https://github.com/magnusdv/forrel/tree/master/data-raw
#'
#' @format Both `FORCE` and `XFORCE` are data frames with the following columns:
#'   * `CHROM`: Chromosome
#'   * `MARKER`: Marker name (rs number)
#'   * `MB`: Physical position in megabases (GRCh38)
#'   * `A1`: First allele
#'   * `A2`: Second allele
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
#' getMap(x, markers = 1:3)
#' getFreqDatabase(x, markers = 1:3)
#'
"FORCE"


#' @rdname FORCE
#' @format NULL
"XFORCE"
