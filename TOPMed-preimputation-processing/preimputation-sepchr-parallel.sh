#!/bin/bash
#SBATCH --time=00:02:00
#SBATCH --account=def-cpolychr
#SBATCH --mem=500MB
#Load required modules
module purge
module load vcftools
module load tabix

#Extract variant data for the current chromosome and save as a new VCF file
vcftools --vcf $PREFIX-FWD_clean.vcf --chr $CHR --recode --recode-INFO-all --out $PREFIX-chr$CHR

#Sort the resulting VCF file and compress it into a gzipped VCF file
vcf-sort $PREFIX-chr$CHR.recode.vcf | bgzip -c > $PREFIX-chr$CHR.vcf.gz
