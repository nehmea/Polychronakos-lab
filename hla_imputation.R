
############################### main parameters ###############################
#prefix for bed files to load. make sure that the .bed, .bim, and .fam files have the smae names
prefix='chr6.dose_filtered-MAF0.01-R0.8_chr628477797-33448354'
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
            'HIBAG_TOPMED-chr6_filtered-MAF0.01-R0.8_chr6-28477797-33448354_rsId.txt',
            sep='\t', col.names = NA)

# plot(density(pred_results[['A']]$value$matching))
# plot.ecdf(pred_results[['A']]$value$matching)
