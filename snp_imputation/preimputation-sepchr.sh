#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=def-cpolychr
#SBATCH -c 4

# Load required modules
module load vcftools
module load tabix

# Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -i|--input) input_file="$2"; shift;;    # Input VCF file path
      -o|--out) output_prefix="$2"; shift;;   # Output file prefix
    esac
    shift
done

# Iterate over chromosomes 1 to 22
# Loop through chromosomes 1 to 22 and perform the following actions for each chromosome:
# Use vcftools to extract variant data for the current chromosome from the input VCF file and save it as a new VCF file with the name $output_prefix-chr$chr.
# Sort the resulting VCF file using vcf-sort.
# Compress the sorted VCF file into a gzipped VCF file with the name $output_prefix-chr$chr.vcf.gz.
for chr in {1..22}
do
	# Extract variant data for the current chromosome and save as a new VCF file
	vcftools --vcf $input_file --chr $chr --recode --recode-INFO-all --out $output_prefix-chr$chr
	
	# Sort the resulting VCF file and compress it into a gzipped VCF file
	vcf-sort $output_prefix-chr$chr.recode.vcf | bgzip -c > $output_prefix-chr$chr.vcf.gz
done


