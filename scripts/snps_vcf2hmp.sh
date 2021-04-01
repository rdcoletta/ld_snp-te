#!/bin/bash
#PBS -l walltime=6:00:00,nodes=1:ppn=1,mem=200gb
#PBS -o /home/hirschc1/della028/projects/ld_snp-te
#PBS -e /home/hirschc1/della028/projects/ld_snp-te
#PBS -V
#PBS -N snps_vcf2hmp
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/ld_snp-te/

# transform to hmp diploid
run_pipeline.pl -Xmx200g -importGuess /home/hirschc1/oconnorc/freebayes_output/SNP_data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf \
                -export analysis/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.hmp.txt \
                -exportType HapmapDiploid
