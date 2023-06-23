
############################### main parameters ###############################
#prefix for bed files to load. make sure that the .bed, .bim, and .fam files have the smae names
prefix='chr6.dose_filtered-MAF0.01-R0.8_chr628477797-33448354_resIds'
# GitHub URL for model params
model_url <- "https://github.com/nehmea/CP-helper-methods/raw/e191710b679443135026d0e8d35b7a8f82553c84/InfiniumGlobal-European-HLA4-hg19.RData"


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
