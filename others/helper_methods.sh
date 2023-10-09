#extract columns and rows from vcf file for specific rsIds
bcftools view --include ID==@grs2_alleles_rsIds.list dbsnp-All-hg38_withChr.vcf.gz --threads 6 | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' -o grs2-alleles-info.txt

#recode plink bed to vcf
plink --bfile [filename prefix] --recode vcf --out [VCF prefix]
gzip --keep [VCF_file]

#vcf to plink bed file
plink -vcf [input_file] --make-bed --out [out_prefix]

#replace ambiguous ids with rsIds in vcf file
bgzip [rsId_pos_file]
tabix -s1 -b2 -e2 "${rsId_file}.gz"
bcftools annotate [vcf_file] --annotations [rsId_pos_file]".gz" -c CHROM,POS,ID -Oz -o [output_file] --threads 6

#extract rows with specific IDs from vcf file
bcftools view -i 'ID=@id_file.txt' dbsnp-All-hg38_withChr.vcf.gz -Oz -o grs2-alleles-info.vcf.gz --threads 6

#extract rows with specific positions from vcf file
#snp_pos_file must include at least CHROM and POS columns
bcftools view [input_file] -R [snp_pos_file] -Oz -o [output_file] --threads 6

#plink score calculation
plink2 --vcf grs2_alleles_rsIds.vcf.gz --score grs2_snp_betas_all_for_plink.txt 1 2 3 header list-variants cols=fid,pheno1,nallele,denom,dosagesum,scoreavgs,scoresums no-mean-imputation ignore-dup-ids --out grs2_snp_score_plink_no-imputation_all-snps

# Concatenate VCF files
bcftools concat [input_files_list] -Oz -o [output_file]

#extract SNPs only
bcftools view --types snps [input_file.vcf.gz] -Oz -o [output_file.vcf.gz]

#recode CHROM field
bcftools annotate --rename-chrs noChr_tochr.txt original.vcf.gz

# download multiple files simultaneously
cat urls.txt | xargs -n 1 -P 2 wget -q