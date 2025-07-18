library(tidyverse)
library(httr2)
library(jsonlite)

xraw = readxl::read_xlsx("data-raw/FORCE-data.xlsx",
                         skip = 2,
                         n_max = 5422,  # skip bad ones at the bottom
                         col_types = "text")

# Ensembl REST API
.getSNPdat = function(rs) {  Sys.sleep(0.2)

  resp = request("https://rest.ensembl.org/variation/human?pops=1") |>
    req_body_json(list(ids = as.list(rs))) |>
    req_headers(Accept = "application/json") |>
    req_perform()

  if (resp_status(resp) == 200)
    resp |> resp_body_string() |> fromJSON()
}

# Download complete data for all snps: Chunk size of 120
n = length(xraw$SNP)
ensData = xraw$SNP |> split(ceiling(seq_len(n) / 120)) |> map(.getSNPdat)

# Any NULLs?
table(sapply(ensData, is.null))
table(lengths(ensData))

# Unlist
ensData2 = ensData |> unname() |> unlist(recursive = F)

# Remove Y
x = xraw |> filter(Chromosome != "Y")

# 1 FORCE SNP missing ("no longer supported"): rs2323964, chr12
setdiff(x$SNP, names(ensData2))

# Set final SNP data
ensData3 = ensData2[intersect(x$SNP, names(ensData2))]

# Save data (for further processing; not for package)
saveRDS(ensData3, "data-raw/ensemblFORCE.rds")

