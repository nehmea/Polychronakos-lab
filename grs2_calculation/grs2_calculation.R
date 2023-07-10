
rm(list=ls())

##############################  parameters ##############################  
#important loci and prob cutoff to filter samples
prob_cutoff = 0.7
important_loci = c('DQA1','DQB1'#, 'A', 'B', 'C', 
                   )


#path to hla imputations. this should be a tab-delimited file. The file should include the following columns:
# locus: HLA locus, make sure the loci names as follows c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")
#sample_id: sample id should be unique for each sample
#allele1, allele2: alleles 1 and 2 for each sample
#prob: probability of imputed HLA haplotype
hla_imputations_file = "../HIBAG_chr6.dose_filtered-MAF0.01-R0.8_chr628477797-33448354_rsId.txt"

#hla and non-hla bed files
hla_bed_prefix='../chr6.dose_filtered-MAF0.01-R0.8_chr628477797-33448354_rsId'
non_hla_prefix='../grs2-non-hla-snps_rsIds'

#beta value files for score calculation
interaction_betas_file = 'https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2-interaction_betas.txt'
drdq_betas_file = 'https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2-dr-dq_betas.txt'
non_drdq_betas_file = 'https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2-non-dr-dq_betas.txt'
non_hla_betas_file = 'https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2-non-hla_betas.txt'

############################### create log file #################################
if(!require(utils)) {
  install.packages(utils)
}

# Load the utils package
library(utils)

# Open the log file for writing
log_file <- "GRS2-log.txt"
sink(log_file, append = TRUE, split=T)

############################### required packages ###############################
# List of required packages
required_packages <- c('ggplot2', 'progress', 'xlsx')

# Check if each required package is installed
missing_packages <- setdiff(required_packages, installed.packages()[, "Package"])

# Install missing packages
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Required Bioconductor packages
required_bioc_packages <- c("HIBAG")

# Check if each required Bioconductor package is installed
missing_bioc_packages <- setdiff(required_bioc_packages, installed.packages()[, "Package"])

# Install missing Bioconductor packages
if (length(missing_bioc_packages) > 0) {
  BiocManager::install(missing_bioc_packages)
}

lapply(c(required_packages, required_bioc_packages), require, character.only = TRUE)

############################### Methods #######################################
#change unwanted characters to avoid errors in reading the file (fam file)
gsub_file = function(file_name, pattern, replacement) {
  fam=readLines(paste0(file_name,".fam"))
  fam=gsub(pattern, replacement, fam)
  writeLines(fam, paste0(file_name,".fam"))
  return(TRUE)
}

#######################  preparing input for GRS2 calculation ####################### 
#reading imputation file
haplotypes = read.table(hla_imputations_file, header=T, sep='\t')

#break if any of the required fields is missing
missing_fields = setdiff(c('locus', 'sample_id', 'allele1', 'allele2', 'prob'), colnames(haplotypes))
if(length(missing_fields > 0)) {
  #cat("\nmissing fields in haplotypes file: ", missing_fields)
  stop(paste0("missing fields in haplotypes file: ", missing_fields))
}

#transforming prob less than prob_cutoff to NA
#included_samples = rownames(hla_imputations_imp_count_prob_cutoff)[hla_imputations_imp_count_prob_cutoff$`count_prob>=0.7` == 4]
haplotypes[haplotypes$prob<prob_cutoff, c('allele1', 'allele2')] <- NA

#perform cutoff to get high-qulaity samples
included_samples = c()
for(sample_id in unique(haplotypes$sample_id)) {
  if(!any(as.vector(haplotypes[haplotypes$sample_id == sample_id & haplotypes$locus %in% important_loci, c('prob')])< prob_cutoff)) {
    included_samples = c(included_samples, sample_id)
  }
}

#haplotypes = hla_imputations[haplotypes$sample_id %in% included_samples, ]


###hla snps ###
#read bed files
hla_genotypes <- hlaBED2Geno(bed.fn=paste0(hla_bed_prefix,".bed"), 
                             fam.fn=paste0(hla_bed_prefix,".fam"),
                             bim.fn=paste0(hla_bed_prefix,".bim"),
                             rm.invalid.allele=TRUE
                             )
#create snp-sample matrix
hla_genotypes_df = as.data.frame(matrix(hla_genotypes$genotype, 
                                        nrow=nrow(hla_genotypes$genotype),
                                        dimnames=list(hla_genotypes$snp.id, hla_genotypes$sample.id)),
                                 check.names=F)
#retain only hq samples
#hla_genotypes_df = hla_genotypes_df[, included_samples]
#delete unwanted object
rm(hla_genotypes)


###non-hla snps###
#remove # chr from fam file to avoid errors
gsub_file(non_hla_prefix, "#", "_")
#read bed files
non_hla_genotypes = hlaBED2Geno(bed.fn=paste0(non_hla_prefix,".bed"), 
                                fam.fn=paste0(non_hla_prefix,".fam"),
                                bim.fn=paste0(non_hla_prefix,".bim"),
                                rm.invalid.allele=TRUE,
                                assembly='hg38',
                                import.chr=""
                                )
#create snp-sample matrix
non_hla_genotypes_df = as.data.frame(matrix(non_hla_genotypes$genotype, 
                                            nrow=nrow(non_hla_genotypes$genotype),
                                            dimnames=list(non_hla_genotypes$snp.id, non_hla_genotypes$sample.id)),
                                     check.names=F)
#retain only hq samples
#non_hla_genotypes_df = non_hla_genotypes_df[, included_samples]
#remove unwanted object
rm(non_hla_genotypes)

###########################  GRS2 calculation ############################ 
#beta value files for score calculation
interaction_betas = data.frame(read.table(interaction_betas_file, header=T))
drdq_betas = read.table(drdq_betas_file, header=T)
non_drdq_betas = read.table(non_drdq_betas_file, header=T)
non_hla_betas = read.table(non_hla_betas_file, header=T)

missing_hla_snps = non_drdq_betas[non_drdq_betas$SNP %in% setdiff(non_drdq_betas$SNP, rownames(hla_genotypes_df)),]
write.table(missing_hla_snps, 'missing_grs2_hla_snps_in_cohort.txt', sep='\t', row.names = F)

missing_non_hla_snps = non_hla_betas[non_hla_betas$SNP %in% setdiff(non_hla_betas$SNP, rownames(non_hla_genotypes_df)),]
write.table(missing_non_hla_snps, 'missing_grs2_non-hla_snps_in_cohort.txt', sep='\t', row.names = F)

grs2_scores = data.frame(matrix(NA, 
                                nrow = length(unique(haplotypes$sample_id)),
                                ncol=4,
                                dimnames=list(unique(haplotypes$sample_id), 
                                              c('interaction_score','drdq_score', 'non_drdq_score', 'non_hla_score'))
                                )
                         )

sample_interaction_betas = interaction_betas
sample_drdq_betas = drdq_betas
sample_non_drdq_betas = non_drdq_betas
sample_non_hla_betas = non_hla_betas

#progress bar
rm(prog_bar)
prog_bar = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                total = length(unique(haplotypes$sample_id)),
                                complete = "=",   # Completion bar character
                                incomplete = "-", # Incomplete bar character
                                current = ">",    # Current bar character
                                clear = FALSE,    # If TRUE, clears the bar when finish
                                width = 100)      # Width of the progress bar

for(sample_id in as.character(unique(haplotypes$sample_id))) {
  
  if(! sample_id %in% included_samples) {
    #grs2_scores[sample_id, 'grs2_score'] = 'low_quality'
    next
    }
  
  #initialize score constituents to zero
  sample_hla_haplotypes = haplotypes[haplotypes$sample_id == sample_id, ]
  
  
  #other hla alleles score
  sample_non_drdq_betas = merge(sample_non_drdq_betas, 
                           data.frame(SNP = sample_non_drdq_betas$SNP,
                                      count = hla_genotypes_df[sample_non_drdq_betas$SNP, as.character(sample_id), drop=F],
                                      check.names = F),
                           by = 'SNP')
  grs2_scores[sample_id, 'non_drdq_score'] = sum(sample_non_drdq_betas$Beta * sample_non_drdq_betas[,sample_id], na.rm = T)
  
  
  #non-hla alleles score
  sample_non_hla_betas = merge(sample_non_hla_betas, 
                               data.frame(SNP = sample_non_hla_betas$SNP,
                                          count = non_hla_genotypes_df[sample_non_hla_betas$SNP, sample_id, drop=F],
                                          check.names = F),
                               by = 'SNP')
  grs2_scores[sample_id, 'non_hla_score'] = sum(sample_non_hla_betas$Beta * sample_non_hla_betas[,sample_id], na.rm = T)
  
  #interaction
  hap1 = sample_hla_haplotypes[sample_hla_haplotypes$locus %in% c('DQA1', 'DQB1'), 'allele1']
  hap2 = sample_hla_haplotypes[sample_hla_haplotypes$locus %in% c('DQA1', 'DQB1'), 'allele2']
  sample_interaction_betas[, sample_id] = (
    c(gsub('X',substr(hap1[1], 5,5),interaction_betas[, c('h1_DQA1')]) %in% hap1[1] &
        interaction_betas[, c('h1_DQB1')] %in% hap1[2] &
        gsub('X',substr(hap2[1], 5,5),interaction_betas[, c('h2_DQA1')]) %in% hap2[1] &
        interaction_betas[, c('h2_DQB1')] %in% hap2[2]
    ) 
    |
      c(gsub('X',substr(hap2[1], 5,5),interaction_betas[, c('h1_DQA1')]) %in% hap2[1] &
          interaction_betas[, c('h1_DQB1')] %in% hap2[2] &
          gsub('X',substr(hap1[1], 5,5),interaction_betas[, c('h2_DQA1')]) %in% hap1[1] &
          interaction_betas[, c('h2_DQB1')] %in% hap1[2]
      )
  )
  
  if(any(sample_interaction_betas[, sample_id])) {
    grs2_scores[sample_id, 'interaction_score'] = interaction_betas$Beta[sample_interaction_betas[, sample_id]]
    grs2_scores[sample_id, 'drdq_score'] = 0
    
    #no interaction
  } else {
    grs2_scores[sample_id, 'interaction_score'] = 0

    sample_drdq_betas_hap1 = as.numeric(gsub('X', substr(hap1[1], 5,5),drdq_betas[, c('DQA1')]) %in% hap1[1] & 
                                          (drdq_betas[, c('DQB1')] %in% c(hap1[2], hap2[2]))) 
    sample_drdq_betas_hap2 = as.numeric(gsub('X',substr(hap2[1], 5,5),drdq_betas[, c('DQA1')]) %in% hap2[1] &
                                          (drdq_betas[, c('DQB1')] %in% c(hap1[2], hap2[2])))
    sample_drdq_betas[,sample_id] = sample_drdq_betas_hap1 + sample_drdq_betas_hap2
    grs2_scores[sample_id, 'drdq_score'] = sum(sample_drdq_betas$Beta * sample_drdq_betas[,sample_id], na.rm = T)
  }
  
  #grs2_scores[sample_id, 'grs2_score'] = sample_interaction_hap_score + sample_drdq_score + sample_hla_score + sample_non_hla_score
  prog_bar$tick()
}

grs2_scores$grs2_score = rowSums(grs2_scores)
#append to EXCEL file
write.xlsx2(grs2_scores, 
            paste0('grs2-scores.xlsx'),
            sheetName = "grs2_scores",
            col.names = TRUE, row.names = T, append = F, showNA=T)


for(score_data in c('sample_interaction_betas', 'sample_drdq_betas', 'sample_non_drdq_betas', 'sample_non_hla_betas')) {
  write.xlsx2(get(score_data), 
              paste0('grs2-scores.xlsx'),
              sheetName = score_data,
              col.names = TRUE, row.names = T, append = T, showNA=T)
}

ggplot(grs2_scores, aes(x='Score', y=grs2_score))+
  geom_violin()+
  geom_boxplot(width=0.1, fill='grey')+ 
  theme_classic()
ggsave('GRS2-score-distribution.png', dpi=300)

percentile_90 = quantile(grs2_scores$grs2_score, probs = 0.9, na.rm = T)
sum(grs2_scores$grs2_score > percentile_90, na.rm=T)

############################## closing log file ################################
# Stop writing to the log file
sink()
sink()

# #TODO: spread haplotypes and add scores and betas to check if they align
# data.frame(sample_id=NA, allel1=NA, allele2=NA,t(interaction_betas[, 1:4]))
# x = haplotypes[haplotypes$locus %in% c('DQA1', 'DQB1'), c('locus','sample_id', 'allele1', 'allele2')]
# rownames(x) = x[, 1]
# t(sample_interaction_betas)
