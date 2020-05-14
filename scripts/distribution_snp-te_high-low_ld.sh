#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=10,mem=130gb
#PBS -o /home/hirschc1/della028/projects/ld_snp-te
#PBS -e /home/hirschc1/della028/projects/ld_snp-te
#PBS -V
#PBS -N distribution_snp-te_high-low_ld_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/ld_snp-te/analysis

module load R

Rscript ~/projects/ld_snp-te/scripts/distribution_snp-te_high-low_ld.R WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.chr${CHR}.ld ld_distribution ${CHR}
