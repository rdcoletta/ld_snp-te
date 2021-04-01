#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("incorrect number of arguments provided.

Usage: Rscript reorder_hmp.R [hmp_file] [header_file]
       ")
}

# assign arguments to variables
hmp_file_snp <- args[1]
hmp_file_te <- args[2]

# hmp_file_snp <- "analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.hmp.txt"
# hmp_file_te <- "analysis/B73_pav_matrix_revision.TA_format.sort.hmp.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### reorder hapmap ----

# load data
hmp_snp <- fread(hmp_file_snp, header = TRUE, data.table = FALSE)
hmp_te <- fread(hmp_file_te, header = TRUE, data.table = FALSE)

# get columns in common
columns_to_keep <- intersect(colnames(hmp_snp), colnames(hmp_te))

# reorder hapmaps
hmp_ordered_snp <- hmp_snp[, columns_to_keep]
hmp_ordered_te <- hmp_te[, columns_to_keep]

# check columns match
cat("SNP and TE hapmaps have the same columns and in same order: ")
all(colnames(hmp_ordered_snp) == colnames(hmp_ordered_te))

# write output
outfile_snp <- gsub("hmp.txt", "reordered.hmp.txt", hmp_file_snp)
outfile_te <- gsub("hmp.txt", "reordered.hmp.txt", hmp_file_te)
fwrite(hmp_ordered_snp, file = outfile_snp, quote = FALSE, row.names = FALSE, sep = "\t", na = NA)
fwrite(hmp_ordered_te, file = outfile_te, quote = FALSE, row.names = FALSE, sep = "\t", na = NA)
