
haplotype_columns= 1:5
jacob_data_original = read.table("T1DGC_3971_GRS2_and_AAb.txt", header = T)
jacob_data_original = jacob_data_original[, haplotype_columns]
colnames(jacob_data_original) = toupper(gsub("hla_", "", colnames(jacob_data_original)))
loci = unique(gsub("_.*","", colnames(jacob_data_original)[2:length(colnames(jacob_data_original))]))


jacob_data = data.frame()
for(locus in loci) {
  locus_alleles = data.frame(jacob_data_original[, c('ANALYTIC_ID')], 
                             locus = locus,
                             jacob_data_original[,grep(paste0("^", locus, "_[0-9]"), colnames(jacob_data_original), value = T)]
                             )
  #locus_alleles = jacob_data_original[, c('ANALYTIC_ID', grep(paste0("^", locus, "_[0-9]"), colnames(jacob_data_original), value = T))]
  colnames(locus_alleles) = c('sample_id', 'locus','allele1', 'allele2')
  jacob_data = rbind(jacob_data, locus_alleles)
}

for(column in c("allele1", "allele2")){
  jacob_data[,column] = sapply(jacob_data[,column], function(x) {
    x= as.character(x)
    if(nchar(x) == 3) { 
      x = paste0("0", x)
      paste0(substr(x, 1, 2), ":", substr(x, 3,4))
    } else {
      paste0(substr(x, 1, 2), ":", substr(x, 3,4))
    }
  })
}


jacob_data$prob = 1

write.table(as.matrix(jacob_data), 'jacob_data_reshaped.txt', sep='\t', row.names = F)


