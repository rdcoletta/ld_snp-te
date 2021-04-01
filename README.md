# LD between SNPs and TEs

Path to the files:

```bash
# TEs
/home/hirschc1/oconnorc/TE_project/TE_calling_WiDiv_panel/B73_TE_calls/SNP_TE_LD_files
#SNPs
/home/hirschc1/oconnorc/freebayes_output/SNP_data
```

Original files were transferred to my MSI account:

```bash
cd ~/projects/ld_snp-te
mkdir {analysis,data,scripts}

# TE info
cp /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE.gff3 data/
cp /home/hirschc1/oconnorc/sarah_te_files/B73.structuralTEv2.2018.12.20.filteredTE_start.in_10bp_rem.nocallTEs.gff3 data/
# TE calls (TEs with many ambiguous calls were already removed)
cp /home/hirschc1/shared/projects/polyTE/B73_pav_matrix_revision.txt data/
```



## TE processing

1. Get middle of TE location and strand, and convert TE calls to hapmap format (absent = `AA`, present = `TT`, ambiguous = `NN`):

> Hapmap format described here: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load

```bash
cut -f 1,4,5,7,9 data/B73.structuralTEv2.2018.12.20.filteredTE.gff3 | tr "=" "\t" | tr ";" "\t" | join -1 6 -2 1 <(sort -k 6,6 -) <(cut -f 9 data/B73.structuralTEv2.2018.12.20.filteredTE_start.in_10bp_rem.nocallTEs.gff3 | sort -k 1,1) | tr " " "\t" | awk '{TEmid=int(($3+$4)/2); print $2"\t"TEmid"\t"$5"\t"$1"\t"$3"\t"$4}' | awk '{if ($3 == ".") print $1"\t"$2"\t+\t"$4"\t"$5"\t"$6"\t"; else print $0}' > analysis/B73_TE_midpoint.strand.txt

# add only TEs that have less than 127 ambiguous calls among the genotypes
head -n 1 data/B73_pav_matrix_revision.txt > analysis/B73_pav_matrix_revision.TA_format.txt
sed 1d data/B73_pav_matrix_revision.txt | sed -e 's/present/TT/g' | sed -e 's/absent/AA/g' | sed -e 's/ambigious/NN/g' >> analysis/B73_pav_matrix_revision.TA_format.txt
```

2. Combine TE data with midpoint/strand file to create a hapmap file:

```bash
head -n 1 analysis/B73_pav_matrix_revision.TA_format.txt | cut --complement -f 1 | awk '{print "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t"$0}' > analysis/B73_pav_matrix_revision.TA_format.hmp.txt

join -1 4 -2 1 <(sort -k 4,4 analysis/B73_TE_midpoint.strand.txt) <(tail -n +2 analysis/B73_pav_matrix_revision.TA_format.txt | sort -k 1,1) | tr " " "\t" | awk '{print $1"\tT\/A\t"$2"\t"$3"\t"$4"\tB73v4\tNA\tNA\tNA\tNA\tNA\t"$0}' | cut -f 1-11,18-527 | grep -v "     B73V4" | sort -k 3,3n -k 4,4n >> analysis/B73_pav_matrix_revision.TA_format.hmp.txt
```

3. Sort newly created hapmap file and transform to hapmap diploid format:

```bash
run_pipeline.pl -SortGenotypeFilePlugin -inputFile analysis/B73_pav_matrix_revision.TA_format.hmp.txt \
                -outputFile analysis/B73_pav_matrix_revision.TA_format.sort.hmp.txt \
                -fileType Hapmap

run_pipeline.pl -importGuess analysis/B73_pav_matrix_revision.TA_format.sort.hmp.txt \
                -export analysis/B73_pav_matrix_revision.TA_format.sort.hmp.txt \
                -exportType HapmapDiploid
```



## SNP Processing

1. Use the SNPs vcf file that was converted to table by Christine:

```bash
# java -jar /home/hirschc1/oconnorc/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantsToTable -V WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf -F CHROM -F POS -GF GT -R ~/maize_refs/B73_chr1-10.fasta -o  WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf.txt

cp /home/hirschc1/oconnorc/freebayes_output/SNP_data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf.txt data/
```

2. Identify SNPs that are within boundaries of TEs:

```bash
tail -n +2 data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf.txt | awk '{print $1"\t"$2"\t"$2}' > data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf.bed

module load bedtools

bedtools intersect -a data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf.bed -b data/B73.structuralTEv2.2018.12.20.filteredTE.gff3 -v -wa | cut -f 1-2 > analysis/WiDiv508_B73v4_allchr_SNPs_notinTEs.txt
# kept 2515350 out of a possible 3146253 sites
```

3. Get vcf file with only SNPs outside of TEs and transform it to hapmap:

```bash
vcftools --vcf /home/hirschc1/oconnorc/freebayes_output/SNP_data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf --out analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr --positions analysis/WiDiv508_B73v4_allchr_SNPs_notinTEs.txt --recode --recode-INFO-all

run_pipeline.pl -Xmx100g -importGuess analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.recode.vcf \
                -export analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.hmp.txt \
                -exportType HapmapDiploid
```



## Calculate LD

1. Merge SNPs and TEs hapmaps:

```bash
# make sure the same genotypes are in both files and in the same order
srun -N 1 --ntasks-per-node=1  --mem-per-cpu=40gb -t 2:00:00 -p interactive --pty bash
module load R/3.6.0
Rscript scripts/reorder_snp-te_hmp.R analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.hmp.txt analysis/B73_pav_matrix_revision.TA_format.sort.hmp.txt

# merge hapmaps
cp analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.reordered.hmp.txt analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt
sed 1d analysis/B73_pav_matrix_revision.TA_format.sort.reordered.hmp.txt >> analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt
```

2. Transform hmp to plink format:

```bash
# note that I needed to sort file before transforming to plink
sbatch scripts/hmp2plk.sh
```

3. Use plink to calculate LD in a 1 Mb window:

```bash
sbatch scripts/plink_te-snp.sh
```

4. Keep R2/D' values only from SNP-TE pairs:

```bash
sbatch scripts/keep_only_te-snp_ld.sh
```



## LD distribution

1. Split data into chromosomes so I can analyze them in parallel later:

```bash
srun -N 1 --ntasks-per-node=10  --mem-per-cpu=100gb -t 2:00:00 -p interactive --pty bash

for chr in {1..10}; do
  zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.ld.gz | head -n 1 > analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.chr$chr.ld
done

# run in parallel
for chr in {1..10}; do
  zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.ld.gz | awk -v chr="$chr" '$1 == chr && $5 == chr' >> analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.chr$chr.ld &
done
```

2. Filter SNP-TE LD files for each chromosome in parallel and generate a file only with high LD and other only only low LD:

```bash
# filter by R2
for chr in {1..10}; do
  sbatch --export=CHR=$chr scripts/distribution_snp-te_high-low_ld_R2.sh
done

# filter by Dprime
for chr in {1..10}; do
  sbatch --export=CHR=$chr scripts/distribution_snp-te_high-low_ld_dprime.sh
done
```

3. Combine different chromosomes into a single file

```bash
# filtered by R2
cp analysis/ld_distribution_R2/plink_results_SNPs-highest-LD-TE_R2_chr1.ld analysis/ld_distribution_R2/plink_results_SNPs-highest-LD-TE_R2.ld
for chr in {2..10}; do
  sed 1d analysis/ld_distribution_R2/plink_results_SNPs-highest-LD-TE_R2_chr${chr}.ld >> analysis/ld_distribution_R2/plink_results_SNPs-highest-LD-TE_R2.ld
done

# filtered by Dprime
cp analysis/ld_distribution_dprime/plink_results_SNPs-highest-LD-TE_Dprime_chr1.ld analysis/ld_distribution_dprime/plink_results_SNPs-highest-LD-TE_Dprime.ld
for chr in {2..10}; do
  sed 1d analysis/ld_distribution_dprime/plink_results_SNPs-highest-LD-TE_Dprime_chr${chr}.ld >> analysis/ld_distribution_dprime/plink_results_SNPs-highest-LD-TE_Dprime.ld
done
```



## Final files

LD between SNPs and TEs only: `analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.ld.gz`
Files with SNPs in highest/lowest LD to TEs: `analysis/ld_distribution_R2` and `analysis/ld_distribution_dprime` folders.
