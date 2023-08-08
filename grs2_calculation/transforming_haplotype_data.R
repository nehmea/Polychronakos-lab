rm(list= ls())

haplotype_columns= c(2:13) #including sample_id column
sample_id_column = 'eID'
add_colon = TRUE
original_data = read.table("Formatted_GRS2_Project_Controls_HLA.txt", header = T)
original_data = original_data[, c(grep(sample_id_column,colnames(original_data)), haplotype_columns)]
colnames(original_data) = toupper(gsub("hla_", "", colnames(original_data)))
loci = unique(gsub("_.*","", colnames(original_data)[2:length(colnames(original_data))]))


reshaped_data = data.frame()
for(locus in loci) {
  locus_alleles = data.frame(original_data[, toupper(sample_id_column)], 
                             locus = locus,
                             original_data[,grep(paste0("^", locus, "_[0-9]"), colnames(original_data), value = T)]
                             )
  #locus_alleles = original_data[, c('ANALYTIC_ID', grep(paste0("^", locus, "_[0-9]"), colnames(original_data), value = T))]
  colnames(locus_alleles) = c('sample_id', 'locus','allele1', 'allele2')
  reshaped_data = rbind(reshaped_data, locus_alleles)
}

if(add_colon){
for(column in c("allele1", "allele2")){
  reshaped_data[,column] = sapply(reshaped_data[,column], function(x) {
    x= as.character(x)
    if(nchar(x) == 3) { 
      x = paste0("0", x)
      paste0(substr(x, 1, 2), ":", substr(x, 3,4))
    } else {
      paste0(substr(x, 1, 2), ":", substr(x, 3,4))
    }
  })
}
}

reshaped_data$prob = 1

write.table(as.matrix(reshaped_data), 'ukbb_haplotypes_for_grs2.txt', sep='\t', row.names = F)

############################## extracting grs2 relevant HLA alleles from haplotypes #######################################
library(stringr)
grs2_hla_alleles = read.table('https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2-non-dr-dq_betas.txt',
                              sep=' ', header=T)

grs2_hla_alleles = grs2_hla_alleles[grepl('.*[*].*',grs2_hla_alleles$Locus),]
grs2_hla_alleles[, c('locus', 'allele')] = str_split_fixed(grs2_hla_alleles$Locus, '[*]', 2)
grs2_hla_alleles$allele = paste0(substr(grs2_hla_alleles$allele, 1, 2), ":", substr(grs2_hla_alleles$allele, 3,4))
grs2_hla_alleles$Locus = paste0(grs2_hla_alleles$locus, '*', grs2_hla_alleles$allele)

allele_data = reshaped_data

library(progress)
rm(prog_bar)
prog_bar = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                            total = length(unique(allele_data$sample_id)),
                            complete = "=",   # Completion bar character
                            incomplete = "-", # Incomplete bar character
                            current = ">",    # Current bar character
                            clear = FALSE,    # If TRUE, clears the bar when finish
                            width = 100)      # Width of the progress bar


for(sample_id in as.character(unique(allele_data$sample_id))) {
  
  sample_alleles= allele_data[allele_data$sample_id == sample_id, ]
  sample_grs2_alleles = data.frame(Locus = grs2_hla_alleles$Locus, 
                                   sample_loci_allele1 = as.numeric(grs2_hla_alleles$Locus %in% paste0(sample_alleles$locus, '*', sample_alleles$allele1)),
                                   sample_loci_allele2 = as.numeric(grs2_hla_alleles$Locus %in% paste0(sample_alleles$locus, '*', sample_alleles$allele2))
                                   )
  
  sample_grs2_alleles$grs2_allele = sample_grs2_alleles$sample_loci_allele1 + sample_grs2_alleles$sample_loci_allele2
  grs2_hla_alleles[, as.character(sample_id)] = as.character(sample_grs2_alleles$grs2_allele)
  prog_bar$tick()
}

write.table(grs2_hla_alleles[, c('SNP', as.character(unique(allele_data$sample_id)))],
            't1dgc_grs2_hla_alleles.txt', sep = '\t', row.names=F)

