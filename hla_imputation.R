
#hlaModelFiles(c("InfiniumGlobal-European-HLA4-hg19.RData"))
# Install the 'httr' package if not already installed
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}

library(httr)

# GitHub URL of the raw RData file
url <- "https://github.com/nehmea/CP-helper-methods/raw/e191710b679443135026d0e8d35b7a8f82553c84/InfiniumGlobal-European-HLA4-hg19.RData"

# Download the RData file to a temporary location
temp_file <- "model_list.RData"
GET(url, write_disk(temp_file))

# Load the RData file
model_list <- get(load(temp_file))

# Clean up the temporary file
unlink(temp_file)

prefix='chr6.dose_filtered-MAF0.01-R0.8_chr628477797-33448354_rsId'
fam=readLines(paste0(prefix,".fam"))
fam=gsub('#', '_', fam)
writeLines(fam, paste0(prefix,".fam"))


library(HIBAG)
genotypes <- hlaBED2Geno(bed.fn=paste0(prefix,".bed"), 
                        fam.fn=paste0(prefix,".fam"),
                        bim.fn=paste0(prefix,".bim"),
                        rm.invalid.allele=TRUE,
                        # assembly='hg38'
                        )
summary(genotypes)


pred_results=list()
for(hla_gene in names(model_list)) {
  print(paste0("HLA-", hla_gene))
  model = hlaModelFromObj(model_list[[hla_gene]])
  pred_guess = hlaPredict(model, genotypes,
                           type="response+prob",
                           match.type = 'RefSNP',
                           same.strand = TRUE
                          )
  print(paste('max_prob =', max(pred_guess$value$prob)))
  pred_results[[hla_gene]]=pred_guess
  rm(pred_guess)
}

saveRDS(pred_results, paste0('HLA-imputations_', prefix,'.rds'))
save.image()
