#!/bin/bash
#PBS -l walltime=10:00:00,nodes=1:ppn=1,mem=110gb
#PBS -o /home/hirschc1/della028/projects/ld_snp-te
#PBS -e /home/hirschc1/della028/projects/ld_snp-te
#PBS -V
#PBS -N plink_te-snp
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/ld_snp-te/

# transform hmp into plink format
plink --file analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.plk --make-founders --r2 gz --ld-window-r2 0 --ld-window 1000000 --ld-window-kb 1000 --out analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined --allow-extra-chr
