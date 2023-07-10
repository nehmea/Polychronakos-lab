#extract columns and rows from vcf file
bcftools view --include ID==@grs2_alleles_rsIds.list ../dbsnp-All-hg38.vcf.gz --threads 6 | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > grs2-alleles-info.txt

#recode plink bed to vcf
plink --bfile [filename prefix] --recode vcf --out [VCF prefix]
gzip --keep [VCF_file]

#replace ambiguous ids with rsIds in vcf file
bgzip [rsId_pos_file]
tabix -p vcf [rsId_pos_file]".gz"
bcftools annotate [vcf_file] --annotations [rsId_pos_file]".gz" -c CHROM,POS,ID -Oz -o [output_file] --threads 6


