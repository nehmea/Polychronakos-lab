rm(list = ls())

##############################  parameters ##############################
# important loci and prob cutoff to filter samples
prob_cutoff <- 0.7
important_loci <- c(
  "DQA1", "DQB1" # , 'A', 'B', 'C',
)

output_grs2_filename = 't1dgc_phased_grs2-scores'
# path to hla imputations. this should be a tab-delimited file. The file should include the following columns:
# locus: HLA locus, make sure the loci names as follows c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")
# sample_id: sample id should be unique for each sample
# allele1, allele2: alleles 1 and 2 for each sample
# prob: probability of imputed HLA haplotype
hla_imputations_file <- "jacob_phased_hap_reshaped.txt"


# hla and non-hla bed files
hla_bed_prefix <- ""
hla_genotypes_df_file <- "t1dgc_grs2_hla_alleles.txt"

non_hla_prefix <- "t1dgc_grs2-snps_3971_rsIds"
non_hla_genotypes_df_file <- ""

# beta value files for score calculation
interaction_betas_file <- "https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2-interaction_betas.txt"
drdq_betas_file <- "https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2-dr-dq_betas.txt"
snp_betas_file <- "https://raw.githubusercontent.com/nehmea/Polychronakos-lab/main/grs2_calculation/grs2_snp_betas_all_for_plink.txt"

############################### create log file #################################
if (!require(utils)) {
  install.packages(utils)
}

# Load the utils package
library(utils)

# Open the log file for writing
log_file <- "GRS2-log.txt"
sink(log_file, append = TRUE, split = T)

############################### required packages ###############################
# List of required packages
required_packages <- c("ggplot2", "progress", "xlsx")

# Check if each required package is installed
missing_packages <- setdiff(required_packages, installed.packages()[, "Package"])

# Install missing packages
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

############################### Methods #######################################

#######################  preparing input for GRS2 calculation #######################
# reading imputation file
haplotypes <- read.table(hla_imputations_file, header = T, sep = "\t")
for (column in c("allele1", "allele2")) {
  haplotypes[, column] <- sapply(haplotypes[, column], function(x) {
    x <- as.character(x)
    if (nchar(x) == 4) {
      x <- paste0("0", x)
    } else {
      x
    }
  })
}

# break if any of the required fields is missing
missing_fields <- setdiff(c("locus", "sample_id", "allele1", "allele2", "prob"), colnames(haplotypes))
if (length(missing_fields > 0)) {
  # cat("\nmissing fields in haplotypes file: ", missing_fields)
  stop(paste0("missing fields in haplotypes file: ", missing_fields))
}

# transforming prob less than prob_cutoff to NA
# included_samples = rownames(hla_imputations_imp_count_prob_cutoff)[hla_imputations_imp_count_prob_cutoff$`count_prob>=0.7` == 4]
haplotypes[haplotypes$prob < prob_cutoff, c("allele1", "allele2")] <- NA

# perform cutoff to get high-qulaity samples
excluded_samples = haplotypes[haplotypes$locus %in% important_loci & (is.na(haplotypes$allele1) | is.na(haplotypes$allele2)),]$sample_id
# haplotypes = hla_imputations[haplotypes$sample_id %in% included_samples, ]

###########################  GRS2 calculation ############################
sample_id_list <- as.character(unique(haplotypes$sample_id))
included_samples = sample_id_list[!sample_id_list %in% excluded_samples]

# beta value files for score calculation
for (filename in c("interaction_betas_file", "drdq_betas_file", "non_drdq_betas_file", "non_hla_betas_file")) {
  filename_value <- get(filename)
  if (filename_value != "" & !is.null(filename_value)) {
    assign(gsub("_file", "", filename), read.table(filename_value, header = T))
  }
}

grs2_scores <- data.frame(matrix(NA,
  nrow = length(sample_id_list),
  ncol = 4,
  dimnames = list(
    sample_id_list,
    c("interaction_score", "drdq_score", "non_drdq_score", "non_hla_score")
  )
))

sample_interaction_betas <- interaction_betas
sample_drdq_betas <- drdq_betas

# progress bar
rm(prog_bar)
prog_bar <- progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = length(unique(haplotypes$sample_id)),
  complete = "=", # Completion bar character
  incomplete = "-", # Incomplete bar character
  current = ">", # Current bar character
  clear = FALSE, # If TRUE, clears the bar when finish
  width = 100
) # Width of the progress bar

for (sample_id in included_samples) {

  # initialize score constituents to zero
  sample_hla_haplotypes <- haplotypes[haplotypes$sample_id == sample_id, ]

  # interaction
  hap1 <- sample_hla_haplotypes[sample_hla_haplotypes$locus %in% c("DQA1", "DQB1"), "allele1"]
  hap2 <- sample_hla_haplotypes[sample_hla_haplotypes$locus %in% c("DQA1", "DQB1"), "allele2"]
  sample_interaction_betas[, sample_id] <- (
    c(gsub("X", substr(hap1[1], 5, 5), interaction_betas[, c("h1_DQA1")]) %in% hap1[1] &
      interaction_betas[, c("h1_DQB1")] %in% hap1[2] &
      gsub("X", substr(hap2[1], 5, 5), interaction_betas[, c("h2_DQA1")]) %in% hap2[1] &
      interaction_betas[, c("h2_DQB1")] %in% hap2[2]) |
      c(gsub("X", substr(hap2[1], 5, 5), interaction_betas[, c("h1_DQA1")]) %in% hap2[1] &
        interaction_betas[, c("h1_DQB1")] %in% hap2[2] &
        gsub("X", substr(hap1[1], 5, 5), interaction_betas[, c("h2_DQA1")]) %in% hap1[1] &
        interaction_betas[, c("h2_DQB1")] %in% hap1[2])
  )

  if (any(sample_interaction_betas[, sample_id])) {
    grs2_scores[sample_id, "interaction_score"] <- interaction_betas$Beta[sample_interaction_betas[, sample_id]]
    grs2_scores[sample_id, "drdq_score"] <- 0

    # no interaction
  } else {
    grs2_scores[sample_id, "interaction_score"] <- 0

    sample_drdq_betas_hap1 <- as.numeric(gsub("X", substr(hap1[1], 5, 5), drdq_betas[, c("DQA1")]) %in% hap1[1] &
      (drdq_betas[, c("DQB1")] %in% c(hap1[2], hap2[2])))
    sample_drdq_betas_hap2 <- as.numeric(gsub("X", substr(hap2[1], 5, 5), drdq_betas[, c("DQA1")]) %in% hap2[1] &
      (drdq_betas[, c("DQB1")] %in% c(hap1[2], hap2[2])))
    sample_drdq_betas[, sample_id] <- sample_drdq_betas_hap1 + sample_drdq_betas_hap2
    grs2_scores[sample_id, "drdq_score"] <- sum(sample_drdq_betas$Beta * sample_drdq_betas[, sample_id], na.rm = T)
  }

  # grs2_scores[sample_id, 'grs2_score'] = sample_interaction_hap_score + sample_drdq_score + sample_hla_score + sample_non_hla_score
  prog_bar$tick()
}


#TODO: calculate grs2_snp_score using plink subprocess
snp_betas

# calculate grs2 score for all samples
grs2_scores$grs2_score <- rowSums(grs2_scores, na.rm = T)
grs2_scores$grs2_score[!complete.cases(grs2_scores)] <- NA

# append to EXCEL file
write.xlsx2(grs2_scores,
  paste0(output_grs2_filename, '.xlsx'),
  sheetName = "grs2_scores",
  col.names = TRUE, row.names = T, append = F, showNA = T
)


for (score_data in c("sample_interaction_betas", "sample_drdq_betas", "sample_non_drdq_betas", "sample_non_hla_betas")) {
  write.xlsx2(get(score_data),
              paste0(output_grs2_filename, '.xlsx'),
    sheetName = score_data,
    col.names = TRUE, row.names = T, append = T, showNA = T
  )
}

library(reshape2)
library(ggplot2)
ggplot(melt(grs2_scores), aes(x = variable, y = value)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "grey") +
  theme_classic()
ggsave(paste0(output_grs2_filename, 'violin_plot.png'), dpi = 300)

#boxplot(grs2_scores)

#percentile_90 <- quantile(grs2_scores$grs2_score, probs = 0.9, na.rm = T)
#sum(grs2_scores$grs2_score > percentile_90, na.rm = T)

############################## closing log file ################################
# Stop writing to the log file
sink()
sink()

# #TODO: spread haplotypes and add scores and betas to check if they align
# data.frame(sample_id=NA, allel1=NA, allele2=NA,t(interaction_betas[, 1:4]))
# x = haplotypes[haplotypes$locus %in% c('DQA1', 'DQB1'), c('locus','sample_id', 'allele1', 'allele2')]
# rownames(x) = x[, 1]
# t(sample_interaction_betas)
