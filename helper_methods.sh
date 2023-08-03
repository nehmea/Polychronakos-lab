#extract columns and rows from vcf file
bcftools view --include ID==@grs2_alleles_rsIds.list ../dbsnp-All-hg38.vcf.gz --threads 6 | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > grs2-alleles-info.txt

#recode plink bed to vcf
plink --bfile [filename prefix] --recode vcf --out [VCF prefix]
gzip --keep [VCF_file]

#replace ambiguous ids with rsIds in vcf file
bgzip [rsId_pos_file]
tabix -p vcf [rsId_pos_file]".gz"
bcftools annotate [vcf_file] --annotations [rsId_pos_file]".gz" -c CHROM,POS,ID -Oz -o [output_file] --threads 6

#extract rows with specific IDs from vcf file
bcftools view -i 'ID=@id_file.txt' [input_file] -Oz -o [output_file]

#plink score calculation
plink2 --vcf t1dgc_grs2-snps_3971_rsIds.vcf.gz --score grs2_snp_betas_all_for_plink.txt 1 2 3 header list-variants cols=fid,pheno1,nallele,denom,dosagesum,scoreavgs,scoresums no-mean-imputation --out GRS2_3971_plink-score_no-imputation_use-alt_all-snps


