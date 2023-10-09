#!/bin/bash
set -u

start=`date +%s`
# Specify the log file name
logfile="logfile_"$start".txt"

# Redirect stdout to both console and log file using tee and unbuffer
exec > >(tee -a "$logfile") 2>&1
exec 2> >(tee -a "$logfile" >&2)
alias echo='echo -e'

# Set default values
Rsq_cutoff=0.8
MAF_cutoff=0.01
perform_phasing=False

# Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
    -r2|--r2_cutoff) Rsq_cutoff="$2"; shift;; 
	-maf|--maf_cutoff) MAF_cutoff="$2"; shift;;
    -vcf|--vcf_file) vcf_file="$2"; shift;;
	-hla_filename|--hla_region_vcf_filename) filtered_vcf="$2"; shift;;
	-p|--phasing) perform_phasing="$2"; shift;;
	-a|--annotation_file) annotation_file="$2"; shift;;
	-gmap|--genomic_map) gmap_file="$2"; shift;;
	-r|--region) region="$2"; shift;;
    esac
    shift
done

# Check if the variable is set
if [ -z "${vcf_file+x}" ]; then
    echo "Error: vcf_file is not defined. please define the initial vcf file using -vcf or --vcf_file"
    exit 1
fi

# Check if the region matches the correct pattern
if [ -z "${region+x}" ]; then
	region=""
	echo "WARNING: region is not assigned, the whole vcf file will be used"
else
	# Regular expression pattern
	pattern="^chr[0-9]+:[0-9]+-[0-9]+$"
	if [[ ! $region =~ $pattern ]]; then
		echo "region is not in the correct format. region should be in the format chr6:28510120-33480577"
		exit 1
	fi
fi

vcf_file=$(basename "$vcf_file" .vcf.gz)
annotation_file=$(basename "$annotation_file" .vcf.gz)
filtered_vcf=$vcf_file"_filtered-MAF"$MAF_cutoff"-R"$Rsq_cutoff"_"$region
vcf_without_chr=$filtered_vcf"_no-chr"
filtered_annotation=$annotation_file"_"$region

if [ -n "${annotation_file}" ]; then
	phasing_input=$filtered_vcf'_rsId'
else
	phasing_input=$filtered_vcf
fi

if [ $perform_phasing == True ]; then
	phased_vcf=$phasing_input"_phased"
fi

# Step 1: Preprocess and QC the imputed genotypes
# Extract HLA region from VCF file and perform QC based on imputation quality and MAF
echo "filtering with MAF >"$MAF_cutoff" AND R2 > "$Rsq_cutoff 

bcftools index $vcf_file".vcf.gz"

if [ -z "${region}" ]; then
	echo "region not defined"
	echo "filtering variants"
	bcftools view $vcf_file".vcf.gz" -i 'R2>'$Rsq_cutoff' && MAF>'$MAF_cutoff  -Oz -o $filtered_vcf".vcf.gz" --threads 12
else
	echo "region = "$region
	echo "filtering and extracting region variants"
	bcftools view $vcf_file".vcf.gz" -r $region -i 'R2>'$Rsq_cutoff' && MAF>'$MAF_cutoff -Oz -o $filtered_vcf".vcf.gz" --threads 12
fi

if [ -z "${annotation_file}"  ]; then
	echo "No annotation file defined; annotation will be skipped."
	
else
	echo "annotating with rsIDs..."
	
	#preparing annotation file
	if [ -z "${region}" ]; then
		echo "region not defined; using whole annotation file"
	else
		echo "extracting region from annotation file"
		region_without_chr=$(echo "$region" | sed 's/^chr//')
		bcftools view -r $region_without_chr $annotation_file".vcf.gz" -Oz -o $filtered_annotation".vcf.gz" --threads 12
	fi
	
	#remove 'chr' from CHROM column in the VCF file
	echo "removing 'chr' from vcf file to be consistent with annotation file"
	zcat $filtered_vcf".vcf.gz" | sed 's/^chr//' | bgzip -c > $vcf_without_chr"-temp.vcf.gz"
	#remove 'chr' from CONTIG_ID in the header
	bcftools view -h $vcf_without_chr"-temp.vcf.gz" | sed 's/^##contig=<ID=chr/##contig=<ID=/' | bcftools reheader -h - $vcf_without_chr"-temp.vcf.gz" -o $vcf_without_chr".vcf.gz"
	rm $vcf_without_chr"-temp.vcf.gz"

	#index annotation and vcf files
	bcftools index $filtered_annotation".vcf.gz"
	bcftools index $vcf_without_chr".vcf.gz"
	
	#perform annotation
	bcftools annotate $vcf_without_chr".vcf.gz" --annotations $filtered_annotation".vcf.gz" -c ID -Oz -o $phasing_input'.vcf.gz'
fi

# The resulting preprocessed and filtered VCF file can then be used for the subsequent steps.

#check if phasing is required
if [ "$perform_phasing" != "True" ]; then
	echo "phasing will not be performed. transforming vcf to bed format..."
    plink --vcf $phasing_input".vcf.gz" --make-bed --out $phasing_input
else
	echo "Performing haplotype phasing..."
    # Perform haplotype phasing using SHAPEIT4.
    bcftools index $phasing_input".vcf.gz" # Indexing vcf file
    shapeit4 --input $phasing_input".vcf.gz" --map $gmap_file --region $region --output $phased_vcf".vcf.gz" --log phased.log --thread 12

    # Transform vcf to bed file, the input format for HIBAG
    plink --vcf $phased_vcf".vcf.gz" --make-bed --out $phased_vcf
fi

#calculate runtime
end=`date +%s`
runtime=`expr $end - $start`
# Convert the runtime to a more readable format
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$((runtime % 60))
printf "Runtime: %02d:%02d:%02d\n" $hours $minutes $seconds

# Exit with a non-zero code to indicate error
exit 1  

