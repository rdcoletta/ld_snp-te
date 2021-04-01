#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100gb
#SBATCH -J hmp2plk
#SBATCH -o /home/hirschc1/della028/projects/ld_snp-te/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/ld_snp-te/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

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
