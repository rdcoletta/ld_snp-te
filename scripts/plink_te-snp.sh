#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=110gb
#SBATCH -J plink_te-snp
#SBATCH -o /home/hirschc1/della028/projects/ld_snp-te/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/ld_snp-te/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/ld_snp-te/

# transform hmp into plink format
plink --file analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window 1000000 --ld-window-kb 1000 --out analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime --allow-extra-chr
