#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=130gb
#SBATCH -J distribution_snp-te_high-low_ld_R2_${CHR}
#SBATCH -o /home/hirschc1/della028/projects/ld_snp-te/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/ld_snp-te/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/ld_snp-te/analysis

module load R

Rscript ~/projects/ld_snp-te/scripts/distribution_snp-te_high-low_ld_R2.R WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.Dprime.SNP-TE-only.chr${CHR}.ld ld_distribution_R2 ${CHR}
