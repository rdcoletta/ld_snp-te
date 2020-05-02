library(data.table)

hmp.file <- "analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.sort.hmp.txt"
header.file <- "analysis/header_hmp_snps-not-in-te.txt"


hmp <- fread(hmp.file, header = TRUE, data.table = FALSE)

header.to.order <- fread(header.file, header = FALSE, data.table = FALSE)
header.to.order <- as.character(header.to.order)

hmp.ordered <- hmp[, header.to.order]

