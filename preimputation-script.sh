#!/bin/bash

# Load required modules
module purge
module load plink/1.9b_6.21-x86_64

# Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      # -i|--input) input_name="$2"; shift;;         # Input file name
      # -o|--output) output_name="$2"; shift;;       # Output file name
      -p|--prefix) PREFIX="$2"; shift;;       # Output file name
	  -str|--strand_file) str_file="$2"; shift;;  # Strand file name
    esac
    shift
done

# Run the plink command to convert the input file to the binary PED format
# filters out all variants and samples with missing call rates exceeding 0.1 (at least 90% call rate)
plink --file $PREFIX --mind 0.1 --geno 0.1 --make-bed --out $PREFIX

# Prepare chromosome, position, and flip files
chr_file=GSAMD-24v3-0-EA_20034606_A1-b37.strand.chr
pos_file=GSAMD-24v3-0-EA_20034606_A1-b37.strand.pos
flip_file=GSAMD-24v3-0-EA_20034606_A1-b37.strand.flip
! test -f $chr_file && cat $str_file | cut -f 1,2 > $chr_file
! test -f $pos_file && cat $str_file | cut -f 1,3 > $pos_file
! test -f $flip_file && cat $str_file | awk '{if ($5=="-") print $0}' | cut -f 1 > $flip_file

# Update chromosome, position, and flip information in the input file
test -f $chr_file && plink --allow-no-sex --bfile $PREFIX --update-map $chr_file --update-chr --make-bed --out TEMP_FILE_1
test -f $pos_file && plink --allow-no-sex --bfile TEMP_FILE_1 --update-map $pos_file --make-bed --out TEMP_FILE_2
test -f $flip_file && plink --allow-no-sex --bfile TEMP_FILE_2 --flip $flip_file --make-bed --out TEMP_FILE_3
test -f $pos_file && plink --allow-no-sex --bfile TEMP_FILE_3 --extract $pos_file --make-bed --out $PREFIX-FWD

# Clean up temporary files
rm -f TEMP_FILE_*

# Create plink text files
plink --bfile $PREFIX-FWD --recode --out $PREFIX-FWD

# Recode plink text files into VCF format
plink --file $PREFIX-FWD --recode vcf --out $PREFIX-FWD

#clean VCF
python cleanVCF.py $PREFIX-FWD


# Loop through chromosomes 1 to 22 
# Submit a separate object for each chromosome that do the following jobs
# Use vcftools to extract variant data for the current chromosome from the input VCF file and save it as a new VCF file with the name $output_prefix-chr$chr.
# Sort the resulting VCF file using vcf-sort.
# Compress the sorted VCF file into a gzipped VCF file with the name $output_prefix-chr$chr.vcf.gz.
for chr in {1..22}
do
	test -f $PREFIX-FWD_clean.vcf && sbatch --export=PREFIX=$PREFIX,CHR=$chr preimputation-sepchr-parallel.sh
done


