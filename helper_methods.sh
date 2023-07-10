#extract columns and rows from vcf file
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' -i 'ID=@non-hla_rsIds.txt' ./dbsnp-All-hg38.vcf.gz > non-hla-allele-position.txt

#recode plink bed to vcf
plink --bfile [filename prefix] --recode vcf --out [VCF prefix]
gzip --keep [VCF_file]

#replace ambiguous ids with rsIds in vcf file
bgzip [rsId_pos_file]
tabix -p vcf [rsId_pos_file]".gz"
bcftools annotate [vcf_file] --annotations [rsId_pos_file]".gz" -c CHROM,POS,ID -Oz -o [output_file]