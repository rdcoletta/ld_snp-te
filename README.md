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
# TE calls
cp /home/hirschc1/oconnorc/TE_project/TE_calling_WiDiv_panel/B73_TE_calls/WiDiv508_ref.B73_3cat.txt data/
```



## TE processing

1. Get middle of TE location and strand, and convert TE calls to hapmap format (absent = `AA`, present = `TT`, ambiguous = `NN`):

> Hapmap format described here: https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load

```bash
cut -f 1,4,5,7,9 data/B73.structuralTEv2.2018.12.20.filteredTE.gff3 | tr "=" "\t" | tr ";" "\t" | join -1 6 -2 1 <(sort -k 6,6 -) <(cut -f 9 data/B73.structuralTEv2.2018.12.20.filteredTE_start.in_10bp_rem.nocallTEs.gff3 | sort -k 1,1) | tr " " "\t" | awk '{TEmid=int(($3+$4)/2); print $2"\t"TEmid"\t"$5"\t"$1"\t"$3"\t"$4}' | awk '{if ($3 == ".") print $1"\t"$2"\t+\t"$4"\t"$5"\t"$6"\t"; else print $0}' > analysis/B73_TE_midpoint.strand.txt

# add only TEs that have less than 127 ambiguous calls among the genotypes
head -n 1 data/WiDiv508_ref.B73_3cat.txt | cut -f 1-509 > analysis/WiDiv508_ref.B73_3cat_low.ambig_TA_call.txt
awk '$511 <= 127' data/WiDiv508_ref.B73_3cat.txt | sed -e 's/Present/TT/g' | sed -e 's/Absent/AA/g' | sed -e 's/ambiguous/NN/g' | cut -f 1-509 >> analysis/WiDiv508_ref.B73_3cat_low.ambig_TA_call.txt
```

2. Combine TE data with midpoint/strand file to create a hapmap file:

```bash
head -n 1 analysis/WiDiv508_ref.B73_3cat_low.ambig_TA_call.txt | cut -f 2-509 | awk '{print "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t"$0}' > analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.hmp.txt

join -1 4 -2 1 <(sort -k 4,4 analysis/B73_TE_midpoint.strand.txt) <(tail -n +2 analysis/WiDiv508_ref.B73_3cat_low.ambig_TA_call.txt | sort -k 1,1) | tr " " "\t" | awk '{print $1"\tT\/A\t"$2"\t"$3"\t"$4"\tB73v4\tNA\tNA\tNA\tNA\tNA\t"$0}' | cut -f 1-11,18-526 | grep -v "     B73V4" | sort -k 3,3n -k 4,4n >> analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.hmp.txt
```

3. Sort newly created hapmap file and transform to hapmap diploid format:

```bash
run_pipeline.pl -SortGenotypeFilePlugin -inputFile analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.hmp.txt \
                -outputFile analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.sort.hmp.txt \
                -fileType Hapmap

run_pipeline.pl -importGuess analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.sort.hmp.txt \
                -export analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.sort.hmp.txt \
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
# kept 2515350 out of a possible 3146253 Sites

module load bedtools

bedtools intersect -a data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf.bed -b data/B73.structuralTEv2.2018.12.20.filteredTE.gff3 -v -wa | cut -f 1-2 > analysis/WiDiv508_B73v4_allchr_SNPs_notinTEs.txt
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
# check if order of genotypes are the same between files
head -n 1 analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.sort.hmp.txt
head -n 1 analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.hmp.txt
# they don't so I have to reorder columns
module load R
Rscript reorder_hmp.R

# merge hapmaps
cp analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.hmp.txt analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt
sed 1d analysis/WiDiv508_ref.B73_3cat_low.ambig_TA.sort.reordered.hmp.txt >> analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt
```

2. Transform hmp to plink format:

```bash
# note that I needed to sort file before transforming to plink
qsub scripts/hmp2plk.sh
```

3. Use plink to calculate LD in a 1 Mb window:

```bash
qsub scripts/plink_te-snp.sh
```

4. Keep R2 values only from SNP-TE pairs:

```bash
qsub scripts/keep_only_te-snp_ld.sh

# compress file for sharing
gzip analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.ld
```



## LD distribution

1. Split data into chromosomes so I can analyze them in parallel later:

```bash
qsub -I -l walltime=2:00:00,nodes=1:ppn=10,mem=100gb

for chr in {1..10}; do
  zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.ld.gz | head -n 1 > analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.chr$chr.ld
done

# run in parallel
for chr in {1..10}; do
  zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.ld.gz | awk -v chr="$chr" '$1 == chr && $4 == chr' >> analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.chr$chr.ld &
done
```

2. Filter SNP-TE LD files for each chromosome in parallel and generate a file only with high LD (R2 > 0.8) and other only only low LD (R2 < 0.2):

```bash
for chr in {1..10}; do
  qsub -v CHR=$chr scripts/distribution_snp-te_high-low_ld.sh
done
```



## Final files

LD between SNPs and TEs only: `analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.ld.gz`
