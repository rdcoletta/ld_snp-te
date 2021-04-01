#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=110gb
#SBATCH -J keep_only_te-snp_ld
#SBATCH -o /home/hirschc1/della028/projects/ld_snp-te/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/ld_snp-te/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/ld_snp-te/

# copy header for filtered file
zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.ld.gz | head -n 1 > analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.ld

# keep only snp and TE r2
zcat analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.ld.gz | awk '$3~/^S[0-9]+_/ && $7~/^[^S][^0-9]/ || $3~/^[^S][^0-9]/ && $7~/^S[0-9]+_/' - >> analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.ld

# compress file for sharing
gzip analysis/WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.ld
