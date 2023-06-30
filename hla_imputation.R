
library(xlsx)

############################### main parameters ###############################
#prefix for bed files to load. make sure that the .bed, .bim, and .fam files have the smae names
prefix='chr6.dose_filtered-MAF0.01-R0.8_chr628477797-33448354_resIds'
# GitHub URL for model params
model_url <- "https://github.com/nehmea/CP-helper-methods/raw/e191710b679443135026d0e8d35b7a8f82553c84/InfiniumGlobal-European-HLA4-hg19.RData"

############################### Methods ###############################
gsub_file = function(file_name, pattern, replacement) {
  fam=readLines(paste0(file_name,".fam"))
  fam=gsub(pattern, replacement, fam)
  writeLines(fam, paste0(file_name,".fam"))
  return(TRUE)
}

############################### loading model params ############################### 
#getting model params available in gitHub repo
# Install the 'httr' package if not already installed
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}

library(httr)

# Download the RData file to a temporary location
temp_file <- "model_list.RData"
GET(model_url, write_disk(temp_file))

# Load the RData file
model_list <- get(load(temp_file))

# Clean up the temporary file
unlink(temp_file)

############################### loading genotypes from bed files ############################### 
##replace '#' with '_' in the .fam file to avoid skipping fields when loading files downstream
fam=readLines(paste0(prefix,".fam"))
fam=gsub('#', '_', fam)
writeLines(fam, paste0(prefix,".fam"))

#loading genotypes
library(HIBAG)
genotypes <- hlaBED2Geno(bed.fn=paste0(prefix,".bed"), 
                        fam.fn=paste0(prefix,".fam"),
                        bim.fn=paste0(prefix,".bim"),
                        rm.invalid.allele=TRUE,
                        # assembly='hg38'
                        )

####################### obtaining HLA predictions using HIBAG ####################### 
#create an empty list to store the results
pred_results=list()

#loop through the different locus models to predict HLA alleles using the HIBAG model
for(hla_gene in names(model_list)) {
  print(paste0("HLA-", hla_gene)) #print current locus being processed
  model = hlaModelFromObj(model_list[[hla_gene]]) #get the HIBAG model
  
  #perform predictions
  pred_guess = hlaPredict(model, genotypes,
                           type="response+prob",
                           match.type = 'RefSNP',
                           same.strand = TRUE
                          )
  
  #print max probability to make sure we are getting good results
  print(paste('max_prob =', max(pred_guess$value$prob)))
  
  #add predictions to results list
  pred_results[[hla_gene]]=pred_guess
  
  #remove current predictions to avoid misinterpretations
  rm(pred_guess)
}

#save prediction results into an rds file
saveRDS(pred_results, paste0('HLA-imputations_', prefix,'.rds'))

#save image for later use
save.image()

#######################  compiling top HLA predictions for each locus ####################### 
#create an empty data frame to store results
best_guess_results=data.frame()

#loop through the each locus, get top guess data, add locus information, and add to final results 
for(loc in names(pred_results)){
  print(loc) #print locus name
  
  # get best guess data and add locus name as new field
  best_guess = data.frame(locus=loc,pred_results[[loc]]$value) 
  
  #add locus predictions to final dataframe
  best_guess_results = rbind(best_guess_results, best_guess)
  
  #delete current data
  rm(best_guess)
}

#write compiled information to txt file
write.table(best_guess_results, 
            paste0('HIBAG_',prefix,'.txt'),
            sep='\t', col.names = NA)

#write to EXCEL file
library(xlsx)
write.xlsx2(best_guess_results, 
            paste0('HIBAG_',prefix,'.xlsx'),
            sheetName = "HIBAG results",
            col.names = TRUE, row.names = F, append = FALSE)
# plot(density(pred_results[['A']]$value$matching))
# plot.ecdf(pred_results[['A']]$value$matching)

########################## statistical summary of results #################################
library(ggplot2)

###statistical summary of posterior probabilities per locus
#create dataframe to store results
statistical_summary=data.frame()
#loop through each locus and append its summary to the results file
for (loc in unique(best_guess_results$locus)) {
  loc_summary=t(as.matrix(summary(best_guess_results[best_guess_results$locus==loc, 'prob'])))
  statistical_summary=rbind(statistical_summary, loc_summary)
}
#assign rownames to the results
rownames(statistical_summary) = unique(best_guess_results$locus)

#append statistical summary to EXCEL file
write.xlsx2(statistical_summary, 
            paste0('HIBAG_',prefix,'.xlsx'),
            sheetName = "prob-summary-per-locus",
            col.names = TRUE, row.names = T, append = T)


####plot posterior probability per locus and save to file
ggplot(best_guess_results, aes(x=locus, y=prob))+
  geom_violin()+ 
  geom_boxplot(width=0.05)+
  labs(title="posterior probability by locus",x="Locus", y = "Posterior Probability")+
  geom_hline(yintercept=0.8, linetype="dashed", color = "red")+
  theme_classic()
ggsave('HIBAG_posterior-probability-by-locus.png', dpi=300)

####get number of samples with post prob below 0.8 or 0.7 for each locus
#prob < 0.8
prob_cutoff_0.8 = as.data.frame.matrix(table(best_guess_results$prob<0.8,
      best_guess_results$locus), row.names = c('prob >= 0.8', 'prob < 0.8'))

#prob < 0.7
prob_cutoff_0.7 = as.data.frame.matrix(table(best_guess_results$prob<0.7,
                                             best_guess_results$locus),
                                       row.names = c('prob >= 0.7', 'prob < 0.7'))

#merge into one dataframe
prob_cutoff_res = rbind(prob_cutoff_0.8, '#','#', prob_cutoff_0.7)

#append to EXCEL file
write.xlsx2(prob_cutoff_res, 
            paste0('HIBAG_',prefix,'.xlsx'),
            sheetName = "prob-cutoff-summary",
            col.names = TRUE, row.names = T, append = T)


####get number of loci with post prob greater than 0.8 or 0.7 for each sample
#get the important locus for T1D score calculation
best_guess_results_imp=best_guess_results[best_guess_results$locus %in% c('A', 'B', 'C', 'DQB1'), ]

#get the number of loci with prob greater than 0.8 or 0.7 for each sample
best_guess_results_imp_count_cutoff_0.8 = t(as.data.frame.matrix(table(best_guess_results_imp$prob<0.8, best_guess_results_imp$sample.id), 
                     row.names = c('count_prob>=0.8', 'count_prob<0.8')))
best_guess_results_imp_count_cutoff_0.7 = t(as.data.frame.matrix(table(best_guess_results_imp$prob<0.7, best_guess_results_imp$sample.id), 
                                                           row.names = c('count_prob>=0.7', 'count_prob<0.7')))
#merge the two dataframes
best_guess_results_imp_count_cutoff = cbind(best_guess_results_imp_count_cutoff_0.8, best_guess_results_imp_count_cutoff_0.7)
best_guess_results_imp_count_cutoff = as.data.frame(best_guess_results_imp_count_cutoff, check.names=F, optional = TRUE)
#append to EXCEL file
write.xlsx2(as.data.frame(best_guess_results_imp_count_cutoff, check.names=F, optional = TRUE), 
            paste0('HIBAG_',prefix,'.xlsx'),
            sheetName = "count_best_guess_imp-locus",
            col.names = TRUE, row.names = T, append = T)

#get the count of samples vs the number of loci with prob greater than 0.8 or 0.7 for each sample
#e.g. what is the number of samples with prob > 0.8 in 1,2,3 or 4 loci
count_samples_prob_threshold = data.frame(
  count_locus=paste0('count_loci=',names(table(factor(best_guess_results_imp_count_cutoff[,1])))),
  count_samples_prob_0.8= as.vector(table(factor(best_guess_results_imp_count_cutoff[,1]))), 
  count_samples_prob_0.7= as.vector(table(factor(best_guess_results_imp_count_cutoff[,3]))),
  row.names=1
  )

#append to EXCEL file
write.xlsx2(as.data.frame(count_samples_prob_threshold, check.names=F, optional = TRUE), 
            paste0('HIBAG_',prefix,'.xlsx'),
            sheetName = "count_samples-vs-prob_threshold",
            col.names = TRUE, row.names = T, append = T)

#######################  preparing input for GRS2 calculation v1 ####################### 
cutoff=0.7
included_samples = rownames(best_guess_results_imp_count_cutoff)[best_guess_results_imp_count_cutoff$`count_prob>=0.7` == 4]
haplotypes = best_guess_results[best_guess_results$sample.id %in% included_samples, ]
haplotypes$locus = paste0('hla_',tolower(haplotypes$locus))
colnames(haplotypes) = gsub('allele', '', colnames(haplotypes))

allele_names <- c()
for (s in unique(haplotypes$locus)) {
  for (n in c(1,2)) {
    combination <- paste(s, n, sep = "_")
    allele_names <- c(allele_names, combination)
  }
}

#######################  preparing input for GRS2 calculation v2 ####################### 

#beta values for score calaculation
interaction_betas = read.table('interaction_betas.txt', header=T)
drdq_betas = read.table('dr-dq_betas.txt', header=T)
non_drdq_betas = read.table('non-dr-dq_betas.txt', header=T)
non_hla_betas = read.table('non-hla_betas.txt', header=T)

#imputed haplotypes
cutoff=0.7
included_samples = rownames(best_guess_results_imp_count_cutoff)[best_guess_results_imp_count_cutoff$`count_prob>=0.7` == 4]
haplotypes = best_guess_results[best_guess_results$sample.id %in% included_samples, ]
haplotypes[haplotypes$prob<cutoff, c('allele1', 'allele2')] <- NA

#hla snp-sample matrix
hla_genotypes = as.data.frame(matrix(genotypes$genotype, 
                                         nrow=nrow(genotypes$genotype),
                                         dimnames=list(genotypes$snp.id, genotypes$sample.id)),
                                  check.names=F)

#non-hla snp-sample matrix
non_hla_prefix='grs2-non-hla-snps_rsIds'
gsub_file(non_hla_prefix, "#", "_")
non_hla_genotypes = hlaBED2Geno(bed.fn=paste0(non_hla_prefix,".bed"), 
                                fam.fn=paste0(non_hla_prefix,".fam"),
                                bim.fn=paste0(non_hla_prefix,".bim"),
                                rm.invalid.allele=TRUE,
                                assembly='hg38',
                                import.chr=""
)
non_hla_genotypes = as.data.frame(matrix(non_hla_genotypes$genotype, 
                            nrow=nrow(non_hla_genotypes$genotype),
                            dimnames=list(non_hla_genotypes$snp.id, non_hla_genotypes$sample.id)),
                            check.names=F)
non_hla_genotypes = non_hla_genotypes[, included_samples]

grs2_scores = data.frame(matrix(NA, nrow = length(unique(haplotypes$sample.id)), dimnames=list(unique(haplotypes$sample.id), 'grs2_score')))
for(sample_name in unique(haplotypes$sample.id)) {
  
  sample_hla_haplotypes = haplotypes[haplotypes$sample.id == sample_name, ]
  sample_interaction_hap_score = 0
  sample_drdq_score = 0
  sample_hla_score = 0
  sample_non_hla_score = 0

  
  
  #other hla alleles
  sample_hla_betas = non_drdq_betas
  sample_hla_betas = merge(sample_hla_betas, 
                                   data.frame(SNP = sample_hla_betas$SNP,
                                              count = hla_genotypes[sample_hla_betas$SNP, sample_name]),
                                   by = 'SNP')
  sample_hla_score = sum(sample_hla_betas$Beta*sample_hla_betas$count, na.rm = T)
  
  #non-hla alleles
  sample_non_hla_betas = non_hla_betas
  sample_non_hla_betas = merge(sample_non_hla_betas, 
                           data.frame(SNP = sample_non_hla_betas$SNP,
                                      count = non_hla_genotypes[sample_non_hla_betas$SNP, sample_name]),
                           by = 'SNP')
  sample_non_hla_score = sum(sample_non_hla_betas$Beta*sample_non_hla_betas$count, na.rm = T)
  
  #interaction
  hap1 = sample_hla_haplotypes[sample_hla_haplotypes$locus %in% c('DQA1', 'DQB1'), 'allele1']
  hap2 = sample_hla_haplotypes[sample_hla_haplotypes$locus %in% c('DQA1', 'DQB1'), 'allele2']
  interaction_haplotype = (
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
  
  if(any(interaction_haplotype)) {
    sample_interaction_hap_score = interaction_betas$Beta[interaction_haplotype]
    sample_drdq_score = 0
    #no interaction
  } else {
    interaction_hap_score = 0
    sample_drdq_betas = drdq_betas
    sample_drdq_betas$hap1 = as.numeric(gsub('X',substr(hap1[1], 5,5),drdq_betas[, c('DQA1')]) %in% hap1[1] & 
         (drdq_betas[, c('DQB1')] %in% hap1[2])) 
    sample_drdq_betas$hap2 = as.numeric(gsub('X',substr(hap2[1], 5,5),drdq_betas[, c('DQA1')]) %in% hap2[1] &
        (drdq_betas[, c('DQB1')] %in% hap2[2]))
    sample_drdq_betas$hap_counts = sample_drdq_betas$hap1 + sample_drdq_betas$hap2
    sample_drdq_score = sum(sample_drdq_betas$Beta * sample_drdq_betas$hap_counts, na.rm = T)
  }
  grs2_scores[sample_name, 'grs2_score'] = sample_interaction_hap_score + sample_drdq_score + sample_hla_score + sample_non_hla_score
}

library(ggplot2)
ggplot(grs2_scores, aes(x='Score',y=grs2_score))+
  geom_violin()+
  geom_boxplot(width=0.1)
