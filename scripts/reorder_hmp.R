#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("incorrect number of arguments provided.

Usage: Rscript reorder_hmp.R [hmp_file] [header_file]
       ")
}

# assign arguments to variables
hmp.file <- args[1]
header.file <- args[2]

# hmp.file <- "analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.sort.hmp.txt"
# header.file <- "analysis/header_hmp_snps-not-in-te.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### reorder hapmap ----

# load data
hmp <- fread(hmp.file, header = TRUE, data.table = FALSE)
header.to.order <- fread(header.file, header = FALSE, data.table = FALSE)
header.to.order <- as.character(header.to.order)

# reorder
hmp.ordered <- hmp[, header.to.order]

# write output
outfile <- gsub("hmp.txt", "reordered.hmp.txt", hmp.file)
fwrite(hmp.ordered, file = outfile, quote = FALSE, row.names = FALSE, sep = "\t", na = NA)
