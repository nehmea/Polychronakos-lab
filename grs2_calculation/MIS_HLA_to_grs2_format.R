library(vcfR)
library(tidyverse)
library(reshape2)

imputed_hla = read.vcfR("../imputed_hla.vcf.gz")
imputed_hla_reshaped = reshape2::melt(extract.gt(imputed_hla, element = "GT")) %>% 
  rename(haplotype = Var1, sample_id = Var2, genotype = value) %>%
  separate(genotype, into=c("allele1", "allele2"), sep="\\|") %>%
  separate(haplotype, into=c("locus", "isoform"), sep="\\*", remove = F) %>%
  #mutate(locus=gsub("\\*.*", "", haplotype)) %>%
  mutate(locus=gsub("HLA_", "", locus)) %>%
  subset(grepl(":", haplotype)) %>%
  subset(allele1 > 0 | allele2 > 0)
  

library(matrixStats)
probs_df = melt(extract.gt(imputed_hla, element = "GP")) %>%
  rename(haplotype = Var1, sample_id = Var2, GP = value) %>%
  separate(GP, into=c("homoREF", "hetero", "homoALT"), sep=",")
probs_df$prob = rowMaxs(apply(as.matrix(probs_df[, c('homoREF', 'hetero', 'homoALT')]), 2, as.numeric))

imputed_hla_reshaped = merge(imputed_hla_reshaped, probs_df, by=c("haplotype", "sample_id"))

allele1 = imputed_hla_reshaped[imputed_hla_reshaped$allele1==1,] %>% 
  mutate(allele1 = isoform) %>%
  select(locus, sample_id, allele1, prob)
allele2 = imputed_hla_reshaped[imputed_hla_reshaped$allele2==1,] %>% 
  mutate(allele2 = isoform)%>%
  select(locus, sample_id, allele2, prob)

haps_df = merge(allele1, allele2, by=c("locus", "sample_id"))
haps_df$prob=rowMeans(haps_df[, c('prob.x', 'prob.y')])

write.table(imputed_hla_reshped, "imputed_hla_reshaped.txt", sep = "\t", row.names = F)
write.table(haps_df, "haplotypes_for_grs2_calculation.txt", sep = "\t", row.names = F)

