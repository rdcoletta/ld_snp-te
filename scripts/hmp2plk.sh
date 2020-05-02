#!/bin/bash
#PBS -l walltime=8:00:00,nodes=1:ppn=1,mem=110gb
#PBS -o /home/hirschc1/della028/projects/ld_snp-te
#PBS -e /home/hirschc1/della028/projects/ld_snp-te
#PBS -V
#PBS -N hmp2plk
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/ld_snp-te/

# sort hmp
run_pipeline.pl -Xmx110g -SortGenotypeFilePlugin -inputFile analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt \
                -outputFile analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt \
                -fileType Hapmap

# transform to hmp diploid
run_pipeline.pl -Xmx110g -importGuess analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt \
                -export analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt \
                -exportType HapmapDiploid

# transform hmp into plink format
run_pipeline.pl -Xmx110g -importGuess analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.hmp.txt \
                -export analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined \
                -exportType Plink
