#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: this script merges SNPs and TEs hapmap files from usda parents to be used in
#              Tassel 5 when projecting TEs into RILs
# 
# Usage: ")
#   quit()
# }

# make sure the correct number of arguments are used
# you should provide 3 arguments
if (length(args) != 3) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript 
       ")
}

# assign arguments to variables
plink.file <- args[1]
out.dir.ld <- args[2]
chr <- args[3]

# setwd("~/projects/ld_snp-te/analysis/")
# plink.file <- "WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.chr1.ld"
# out.dir.ld <- "ld_distribution"
# chr <- 10

#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParalell")


if (detectCores() > 10) {
  num.cores <- 10
} else {
  num.cores <- detectCores()
}

if (!dir.exists(out.dir.ld)) dir.create(out.dir.ld)

#### subsample ----

# load one chr at a time
plink.file.chr <- gsub("chr[0-9]+", paste0("chr", chr), plink.file, perl = TRUE)
# open table with LD among markers
LD_results <- fread(plink.file.chr, header = TRUE, data.table = FALSE)

# get TE names
TEs.chr <- append(LD_results[grep("^S[0-9]+_", LD_results[, "SNP_A"], perl = TRUE, invert = TRUE), "SNP_A"],
                  LD_results[grep("^S[0-9]+_", LD_results[, "SNP_B"], perl = TRUE, invert = TRUE), "SNP_B"])
TEs.chr <- unique(TEs.chr)

# add column with distance between te and snp
LD_results$dist_to_te <- LD_results[, 5] - LD_results[, 2]


cat("subsetting only SNPs with highest LD to TE\n")


# create empty dataset
LD_results_highest <- data.frame(matrix(nrow = 0, ncol = NCOL(LD_results)), stringsAsFactors = FALSE)
colnames(LD_results_highest) <- colnames(LD_results)

# # create vector to make sure the same snp doesn't get picked twice
# snps.already.in.ld <- c()

# get closest (highest LD) snps
for (te in TEs.chr) {
  
  # subset LD results to have only the TE being parsed
  snps_LD_with_te <- LD_results[which(LD_results[, "SNP_A"] == te | LD_results[, "SNP_B"] == te), ]

  # if (length(snps.already.in.ld) > 0) {
  #   snps_LD_with_te <- snps_LD_with_te[which(!snps_LD_with_te[, "SNP_A"] %in% snps.already.in.ld & !snps_LD_with_te[, "SNP_B"] %in% snps.already.in.ld), ]
  # }
  
  
  if (NROW(snps_LD_with_te) > 0) {
    
    # select only SNP with highest LD with that TE
    snps_LD_with_te <- snps_LD_with_te[which(snps_LD_with_te[, "R2"] == max(snps_LD_with_te[, "R2"])), ]
    # if there are more than one SNP with the same R2, get the closest one to the TE
    if (NROW(snps_LD_with_te) > 1) {
      snps_LD_with_te <- snps_LD_with_te[which(snps_LD_with_te[, "dist_to_te"] == min(snps_LD_with_te[, "dist_to_te"])), ]
    }
    
    ## get closest snp
    # snp.selected <- apply(snps_LD_with_te, MARGIN = 1, function(row) {
    #   marker1 <- row["SNP_A"]
    #   marker2 <- row["SNP_B"]
    #   if (grepl(paste0("^S", chr, "_"), marker1)) {
    #     return(marker1)
    #   } else {
    #     return(marker2)
    #   }
    # })
    # snps.already.in.ld <- append(snps.already.in.ld, as.character(snp.selected))
    
    # add closest SNP in LD with TE into new df
    LD_results_highest <- rbind(LD_results_highest, snps_LD_with_te)
    
  }
}

# remove duplicates
LD_results_highest <- LD_results_highest[!duplicated(LD_results_highest[, c("SNP_A", "SNP_B")]), ]

# write filtered ld table
outfile.highest <- paste0(out.dir.ld, "/plink_results_SNPs-highest-LD-TE_chr", chr, ".ld")
fwrite(LD_results_highest, outfile.highest, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)



cat("subsetting only SNPs not in LD to TE\n")

# first get names of SNPs that are in LD with an TE (R2>0.2)
snps.in.ld <- mclapply(1:NROW(LD_results), function(row, LD_results) {
  
  marker1 <- LD_results[row, "SNP_A"]
  marker2 <- LD_results[row, "SNP_B"]
  r2 <- as.numeric(LD_results[row, "R2"])
  
  if (grepl(paste0("^S", chr, "_"), marker1)) {
    snp <- marker1
  } else {
    snp <- marker2
  }
  
  if (r2 >= 0.2) return(snp)
  
}, LD_results, mc.cores = num.cores)

snps.in.ld.vector <- do.call(c, snps.in.ld)
snps.in.ld.vector <- snps.in.ld.vector[!duplicated(snps.in.ld.vector)]

# exclude such SNPs
snps.to.exclude <- which(!LD_results[, "SNP_A"] %in% snps.in.ld.vector & !LD_results[, "SNP_B"] %in% snps.in.ld.vector)
LD_results_lowest <- LD_results[snps.to.exclude, ]

# remove duplicates
LD_results_lowest <- LD_results_lowest[!duplicated(LD_results_lowest[, c("SNP_A", "SNP_B")]), ]

# write filtered ld table
outfile.lowest <- paste0(out.dir.ld, "/plink_results_SNPs-lowest-LD-TE_chr", chr, ".ld")
fwrite(LD_results_lowest, outfile.lowest, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)
