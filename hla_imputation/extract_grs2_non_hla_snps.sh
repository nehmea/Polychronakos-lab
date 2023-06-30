#!/bin/bash

max_cores=2
remove_files=false  # Default value for the --remove flag
# Directory paths
vcf_zip_dir="/mnt/d/WORK/Polychronakos/03-Work/02-JDRF/Final_All_samples_JDRF_DM-P01/PLINK_220323_0442/callrate-0.98/TOPMED-imputations/imputations"
output_dir="output_vcf"

# filenames
snp_file="non-hla-alleles-region-file.txt"
bed_output="concatenated"
#rsId_file="non-hla-alleles-annotation-file.txt"

# List of chromosomes
chromosomes=("1" "2" "4" "6" "7" "9" "10" "11" "12" "13" "14" "15" "16" "18" "19" "20" "21" "22")
#chromosomes=(seq 1 22)
#chromosomes=("22")

# Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
	-chr|chromosomes_list) chromosomes="$2"; shift;; 
    -r|--remove) remove_files=true; shift;; 
	-zipDir|--vcfDir) vcf_zip_dir="$2"; shift;;
    -out|--output_dir) output_dir="$2"; shift;;
	-snp|--snp_pos) snp_file="$2"; shift;;
	-rsid|--rsid_file) rsId_file="$2"; shift;;
	-bed|--output_bed) bed_output="$2"; shift;;
	-cores|--workers) max_cores="$2"; shift;;
    esac
    shift
done

# Create the output directory if it does not exist
sudo mkdir "./${output_dir}"
sudo mkdir "./${output_dir}/bed_files"




# Loop through each chromosome
process_chromosome() {
	local chr_num="$1"
	local vcf_zip_dir="$2"
	local snp_file="$3"
	local output_dir="$4"

	echo "processing chromsome ${chr_num}"
		
	# Zip file
    zip_file="${vcf_zip_dir}/chr_${chr_num}.zip"

    # Extract VCF file from the zip
    vcf_gz_file="chr${chr_num}.dose.vcf.gz"
    unzip -n -P "5cK9AdADTJluo7" "${zip_file}" "chr${chr_num}.dose.vcf.gz" -d "${output_dir}"

    #Check if the extraction was successful
    if [ -f "${output_dir}/chr${chr_num}.dose.vcf.gz" ]; then
		
		#create index files
		bcftools index "${output_dir}/chr${chr_num}.dose.vcf.gz"
		
        # Extract the snps with the specified regions in snp_file
        bcftools view "${output_dir}/chr${chr_num}.dose.vcf.gz" -R "${snp_file}" -Oz -o "${output_dir}/chr${chr_num}_filtered.vcf.gz"
		
		#create bed files
		#plink --vcf "${output_dir}/chr${chr_num}_filtered.vcf.gz" --make-bed --out "${output_dir}/bed_files/chr${chr_num}_filtered"
		
        # Delete the .vcf.gz file if the --remove flag is present
        if [ "$remove_files" = true ]; then
            rm "${output_dir}/chr${chr_num}.dose.vcf.gz"
        fi

        echo "Chromosome ${chr_num} processed successfully."
    else
		echo "${output_dir}/chr${chr_num}.dose.vcf.gz"
        echo "Failed to extract VCF file for chromosome ${chr_num}."
    fi
}

# Export the function to make it accessible to parallel
export -f process_chromosome

parallel --jobs "$max_cores" process_chromosome ::: "${chromosomes[@]}" ::: "$vcf_zip_dir" ::: "$snp_file" ::: "$output_dir"

# Wait for all background processes to finish
wait

# List all VCF files in the directory
vcf_files=$(ls "${output_dir}"/*_filtered.vcf.gz)

# Concatenate VCF files
bcftools concat ${vcf_files} -Oz -o "${output_dir}/${bed_output}.vcf.gz"

if [ -z "${rsId_file+x}" ]; then
	#create bed files
	plink --vcf "${output_dir}/${bed_output}.vcf.gz" --make-bed --out "${output_dir}/${bed_output}"
else
	echo "WARNING: annotating with rsIDs"
	#perform annotation
	bgzip "$rsId_file"
	tabix -p vcf "${rsId_file}.gz"
	bcftools annotate "${output_dir}/${bed_output}.vcf.gz" --annotations "${rsId_file}.gz" -c CHROM,POS,ID -Oz -o "${output_dir}/${bed_output}_rsIds.vcf.gz"
		
	#create bed files
	plink --vcf "${output_dir}/${bed_output}_rsIds.vcf.gz" --make-bed --out "${output_dir}/${bed_output}_rsIds"
fi

