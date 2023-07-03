
rm(list=ls())
############################### required packages ###############################
# List of required packages
required_packages <- c("xlsx", "httr", 'ggplot2')

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

############################### main parameters ###############################
#prefix for bed files to load. make sure that the .bed, .bim, and .fam files have the same names
hla_bed_prefix='chr6.dose_filtered-MAF0.01-R0.8_chr628477797-33448354_rsId'

#important loci
important_loci = c('A', 'B', 'C', 'DQA1','DQB1')

# GitHub URL for model params
model_url <- "https://github.com/nehmea/CP-helper-methods/raw/e191710b679443135026d0e8d35b7a8f82553c84/InfiniumGlobal-European-HLA4-hg19.RData"


############################### Methods #######################################
#change unwanted characters to avoid errors in reading the file (fam file)
gsub_file = function(file_name, pattern, replacement) {
  fam=readLines(paste0(file_name,".fam"))
  fam=gsub(pattern, replacement, fam)
  writeLines(fam, paste0(file_name,".fam"))
  return(TRUE)
}

impute_hla = function(model_list, hla_genotypes) {
  
  require(HIBAG)
  
  #create an empty list to store the results
  pred_results=list()
  
  #loop through the different locus models to predict HLA alleles using the HIBAG model
  for(hla_gene in names(model_list)) {
    cat(paste0("\nHLA-", hla_gene, "\n")) #print current locus being processed
    model = hlaModelFromObj(model_list[[hla_gene]]) #get the HIBAG model
    
    #perform predictions
    pred_guess = hlaPredict(model, hla_genotypes,
                            type="response+prob",
                            match.type = 'RefSNP',
                            same.strand = TRUE
                            )
    
    #print max probability to make sure we are getting good results
    cat(paste0('\nmax_prob = ', max(pred_guess$value$prob), "\n"))
    
    #add predictions to results list
    pred_results[[hla_gene]]=pred_guess
    
    #remove current predictions to avoid misinterpretations
    rm(pred_guess)
    
  }
  
  #return results
  return(pred_results)
}

############################### loading model params ############################### 
cat("\n########## loading model parameters ##########\n")

#getting model params available in gitHub repo
# Download the RData file to a temporary location
temp_file <- "model_list.RData"
GET(model_url, write_disk(temp_file))

# Load the RData file
model_list <- get(load(temp_file))

# Clean up the temporary file
unlink(temp_file)

############################### loading genotypes from bed files ############################### 
cat("\n########## loading hla genotypes ##########\n")
##replace '#' with '_' in the .fam file to avoid skipping fields when loading files downstream
gsub_file(hla_bed_prefix, "#", "_")

#loading genotypes

hla_genotypes <- hlaBED2Geno(bed.fn=paste0(hla_bed_prefix,".bed"), 
                        fam.fn=paste0(hla_bed_prefix,".fam"),
                        bim.fn=paste0(hla_bed_prefix,".bim"),
                        rm.invalid.allele=TRUE,
                        # assembly='hg38'
                        )

####################### obtaining HLA predictions using HIBAG ####################### 
cat("\n########## HLA imputation ##########\n")
#create an empty list to store the results
pred_results=list()

#loop through the different locus models to predict HLA alleles using the HIBAG model
pred_results = impute_hla(model_list, hla_genotypes)

#save prediction results into an rds file
saveRDS(pred_results, paste0('HIBAG-imputations_', hla_bed_prefix,'.rds'))

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
colnames(best_guess_results) = gsub("[.]", "_", colnames(best_guess_results))
write.table(best_guess_results, 
            paste0('HIBAG_', hla_bed_prefix,'.txt'),
            sep='\t', col.names = NA)

#write to EXCEL file
library(xlsx)
write.xlsx2(best_guess_results, 
            paste0('HIBAG_', hla_bed_prefix,'.xlsx'),
            sheetName = "HIBAG results",
            col.names = TRUE, row.names = F, append = FALSE)
# plot(density(pred_results[['A']]$value$matching))
# plot.ecdf(pred_results[['A']]$value$matching)

########################## statistical summary of results #################################
cat("\n########## Summarizing imputation results ##########\n")
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
            paste0('HIBAG_', hla_bed_prefix,'.xlsx'),
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
            paste0('HIBAG_', hla_bed_prefix,'.xlsx'),
            sheetName = "prob-cutoff-summary",
            col.names = TRUE, row.names = T, append = T)


####get number of loci with post prob greater than 0.8 or 0.7 for each sample
#get the important locus for T1D score calculation
best_guess_results_imp=best_guess_results[best_guess_results$locus %in% important_loci, ]

#get the number of loci with prob greater than 0.8 or 0.7 for each sample
best_guess_results_imp_count_cutoff_0.8 = t(as.data.frame.matrix(table(best_guess_results_imp$prob<0.8, best_guess_results_imp$sample_id), 
                     row.names = c('count_prob>=0.8', 'count_prob<0.8')))
best_guess_results_imp_count_cutoff_0.7 = t(as.data.frame.matrix(table(best_guess_results_imp$prob<0.7, best_guess_results_imp$sample_id), 
                                                           row.names = c('count_prob>=0.7', 'count_prob<0.7')))
#merge the two dataframes
best_guess_results_imp_count_cutoff = cbind(best_guess_results_imp_count_cutoff_0.8, best_guess_results_imp_count_cutoff_0.7)
best_guess_results_imp_count_cutoff = as.data.frame(best_guess_results_imp_count_cutoff, check.names=F, optional = TRUE)

#append to EXCEL file
write.xlsx2(as.data.frame(best_guess_results_imp_count_cutoff, check.names=F, optional = TRUE), 
            paste0('HIBAG_', hla_bed_prefix,'.xlsx'),
            sheetName = "count-important-loci-above-threshold-per-sample",
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
            paste0('HIBAG_', hla_bed_prefix,'.xlsx'),
            sheetName = "count-samples-with-high-qc-for-imp-loci",
            col.names = TRUE, row.names = T, append = T)
