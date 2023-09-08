#!/bin/bash

max_cores=4
remove_files=false  # Default value for the --remove flag
password=""
# Directory paths
vcf_zip_dir=""
vcf_dir="."

# filenames
snp_file="grs2-alleles-info-hg19.txt"
bed_output="concatenated"
rsId_file="grs2-alleles-info-hg19.txt"

# List of chromosomes
#chromosomes=("1" "2" "4" "6" "7" "9" "10" "11" "12" "13" "14" "15" "16" "18" "19" "20" "21" "22")
chromosomes=$(seq 1 22)
#chromosomes=("22")

# Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
	-chr|chromosomes_list) chromosomes="$2"; shift;; 
    -r|--remove) remove_files=true; shift;; 
	-vcf|--vcfDir) vcf_dir="$2"; shift;;
	-zip|--zip_dir) vcf_zip_dir="$2"; shift;;
	-snp|--region_file) snp_file="$2"; shift;;
	-rsid|--rsid_file) rsId_file="$2"; shift;;
	-bed|--output_bed) bed_output="$2"; shift;;
	-cores|--workers) max_cores="$2"; shift;;
	-p|--password) password="$2"; shift;;
    esac
    shift
done

# Create the output directory if it does not exist
# check if vcf_zip_dir provided
if [ -n "$vcf_zip_dir" ] && [ -n "${vcf_zip_dir// }" ]; then
	vcf_dir="./vcf_dir"
	sudo mkdir "${vcf_dir}"
	sudo mkdir "${vcf_dir}/filtered"
	
#if not, check if vcf_dir provided
else
	#make directory for filtered output
	if [ "$vcf_dir" == "." ]; then
		echo "vcf directory set to current directory."
	fi
	
	ls_result=$(ls $vcf_dir/*.vcf.gz)
	if [ -n "$ls_result" ]; then
		sudo mkdir "${vcf_dir}/filtered"
	else
		echo "Error: no vcf file found in vcf directory."
		exit 1
	fi
fi

if [ -f "${snp_file}" ]; then
	chromosomes=$(awk -F'\t' 'NR==1 {next} !seen[$1]++ {print $1}' $snp_file)
else
	echo "Error: SNPs file ${snp_file} does not exist."
	exit 1
fi
	
# Loop through each chromosome
process_chromosome() {
	local chr_num_initial="$1"
	local vcf_zip_dir="$2"
	local snp_file="$3"
	local vcf_dir="$4"
	local max_cores="$5"
	
	chr_num="${chr_num_initial#chr}"

	echo "processing chromsome ${chr_num}"
		
	# Zip file
	if [ -n "$vcf_zip_dir" ] && [ -n "${vcf_zip_dir// }" ]; then
		zip_file=$(ls $vcf_zip_dir/*chr_$chr_num*.zip)

		# Extract VCF file from the zip
		if [ -e "$zip_file" ]; then
			#vcf_gz_file="*chr${chr_num}*.vcf.gz"
			unzip -n -P $password -d "${vcf_dir}"
		fi
	fi
	
	vcf_gz_file=$(ls $vcf_dir/*chr${chr_num}[!0-9]*.vcf.gz)
	#Check if vcf.gz exists
	if [ -n "$vcf_gz_file" ]; then
		echo "vcf_gz_file = $vcf_gz_file"
		
		#create index files
		# Check if both .tbi and .csi index files do not exist
		#if [ ! -e "${vcf_dir}/${vcf_gz_file}.tbi" ] && [ ! -e "${vcf_dir}/${vcf_gz_file}.csi" ]; then
		bcftools index $vcf_gz_file --threads $max_cores
		#fi
		
		# Extract the snps with the specified regions in snp_file
		bcftools view $vcf_gz_file -R "${snp_file}" -Oz -o "${vcf_dir}/filtered/chr${chr_num}_filtered.vcf.gz"
		
		# Delete the .vcf.gz file if the --remove flag is present
		if [ "$remove_files" = true ] && [ -n "$vcf_zip_dir" ] && [ -n "${vcf_zip_dir// }" ]; then
			rm $vcf_dir/*chr$chr_num*.gz
		fi

		echo "Chromosome ${chr_num} processed successfully."
	else
		echo "$vcf_gz_file not found or failed to extract from zip"
	fi
}

# Export the function to make it accessible to parallel
export -f process_chromosome

parallel --jobs "$max_cores" process_chromosome ::: "${chromosomes[@]}" ::: "$vcf_zip_dir" ::: "$snp_file" ::: "$vcf_dir" ::: "$max_cores"

# Wait for all background processes to finish
wait

# Concatenate VCF files
# List all VCF files in the directory
vcf_files=$(ls $vcf_dir/filtered/*_filtered.vcf.gz)

#check if files exist
if [ -n "$vcf_files" ]; then

	#concatenate
	bcftools concat ${vcf_files} -Oz -o "${vcf_dir}/filtered/concatenated_snps.vcf.gz"
	
	#annotating with rsIds
	if [ -z "${rsId_file+x}" ]; then
		echo "No annotation file provided."
	else
		echo "WARNING: annotating with rsIDs"
		#perform annotation
		bgzip "$rsId_file"
		tabix -s1 -b2 -e2 "${rsId_file}.gz"
		bcftools annotate "${vcf_dir}/filtered/concatenated_snps.vcf.gz" --annotations "${rsId_file}.gz" -c CHROM,POS,ID -Oz -o "${vcf_dir}/filtered/concatenated_snps_rsIds.vcf.gz"
	fi
fi

#remove indels
bcftools view --types snps "${vcf_dir}/filtered/concatenated_snps_rsIds.vcf.gz" -Oz -o "${vcf_dir}/filtered/concatenated_snps_rsIds_snvs.vcf.gz"
