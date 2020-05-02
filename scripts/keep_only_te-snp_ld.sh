#!/bin/bash
#PBS -l walltime=10:00:00,nodes=1:ppn=1,mem=110gb
#PBS -o /home/hirschc1/della028/projects/ld_snp-te
#PBS -e /home/hirschc1/della028/projects/ld_snp-te
#PBS -V
#PBS -N keep_only_te-snp_ld
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/ld_snp-te/

# copy header for filtered file
zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.ld.gz | head -n 1 > analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.ld

# keep only snp and sv r2 (excluding translocations)
zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.ld.gz | awk '$3~/^S[0-9]+_/ && $6~/^[^S][^0-9]/ || $3~/^[^S][^0-9]/ && $6~/^S[0-9]+_/' - >> analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.ld
