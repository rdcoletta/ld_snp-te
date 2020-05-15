#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: this script merges SNPs and SVs hapmap files from usda parents to be used in
#              Tassel 5 when projecting SVs into RILs
#
# Usage: ")
#   quit()
# }

# make sure the correct number of arguments are used
# you should provide 2 arguments
if (length(args) != 2) {
  stop("incorrect number of arguments provided.
       Usage:
       ")
}

# assign arguments to variables
ld.folder <- args[1]
plot_name <- args[2]

# setwd("~/projects/ld_snp-te/analysis/")
# ld.folder = "ld_distribution"
# plot_name <- "dist-LD_SNPs-TEs_high.png"
# plot_name <- "dist-LD_SNPs-TEs_low.png"




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### plot ----

# get filenames
ld.files <- list.files(path = ld.folder, pattern = ".ld", full.names = TRUE)

# select only filenames from high or low LD
if (grepl("high", plot_name)) ld.files <- ld.files[grepl("high", ld.files)]
if (grepl("low", plot_name)) ld.files <- ld.files[grepl("low", ld.files)]

# merge info from chromosomes
ld.df <- data.frame(stringsAsFactors = FALSE)
for (chr in 1:10) {
  ld.chr <- ld.files[grep(paste0("chr", chr, "."), ld.files, fixed = TRUE)]
  ld.chr <- fread(ld.chr, header = TRUE, data.table = FALSE)
  ld.df <- rbind(ld.df, ld.chr)
}

# get the plot subtitle depending on which file i'm analyzing
if (grepl("high", plot_name)) plot_subtitle <- "Subset of SNPs in high LD to TEs"
if (grepl("low", plot_name)) plot_subtitle <- "Subset of SNPs in low LD to TEs"


# distribution of r2 of SNPs in LD with SVs
dist.ld <- ggplot(ld.df, aes(x = R2)) +
  geom_histogram(fill = "#900721", binwidth = 0.005) +
  labs(x = bquote("LD"~(r^2)),
       y = "Count",
       subtitle = plot_subtitle) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

# print how many snps are (or are not) in LD with TE
if (grepl("high", plot_name)) cat(NROW(ld.df) - 1, " closest SNPs in LD with a TE\n", sep = "")
if (grepl("low", plot_name)) cat(NROW(ld.df) - 1, " SNPs not in LD with a TE\n", sep = "")

ggsave(plot = dist.ld, filename = paste0(ld.folder, "/", plot_name), device = "png")
