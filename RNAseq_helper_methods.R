#The mergeCounts function in R is designed to merge counts data from multiple files into a single output. 
#It reads the gene and counts columns from each file and combines them based on matching gene column values.

# Parameters
    # files: A vector of file paths representing the input files containing the counts data.
    # header: A boolean indicating whether the input files have a header row. Default is TRUE.
    # geneCol: The column index or name specifying the gene column in the input files.
    # countsCol: The column index or name specifying the counts column in the input files.
    # sampleNames: An optional vector of sample names to be used as column names in the merged output. 
      #If not provided or the length is not equal to the number of columns in the output, default column names will be used.
    # Output
    # The function returns a merged counts matrix or data frame where each row represents a gene and each column represents a sample.
    #The row names are set as the gene column values,
    #and the column names are determined based on the sampleNames argument or using a default naming scheme.
  # 
  # Progress Bar
  # The function utilizes the progress package to display a progress bar during the merging process. 
    #The progress bar provides information about the completion percentage, elapsed time, and estimated time remaining.
# 
# Usage Example
    # Copy code
    # # Define file paths and column indices/names
    # files <- c("file1.txt", "file2.txt", "file3.txt")
    # geneCol <- "Gene"
    # countsCol <- "Counts"
    # 
    # # Merge counts data
    # mergedData <- mergeCounts(files, geneCol = geneCol, countsCol = countsCol)
# 
# # Display the merged data
# print(mergedData)
# In the above example, mergeCounts is called with three file paths, specifying the gene and counts columns. 
#The resulting merged counts data is stored in the mergedData variable and then printed to the console.
# 
# Note: Ensure that the progress package is installed and loaded before using the mergeCounts function.
mergeCounts = function(files, header=T, geneCol, countsCol, sampleNames=NULL){
  require(progress)
  output = read.table(files[1], header = header)[,c(geneCol,countsCol), drop=F]
  progress_bar = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                  total = length(files),
                                  complete = "=",   # Completion bar character
                                  incomplete = "-", # Incomplete bar character
                                  current = ">",    # Current bar character
                                  clear = FALSE,    # If TRUE, clears the bar when finish
                                  width = 100)      # Width of the progress bar
  for(file in files[2:length(files)]) {
    newData = read.table(file, header = header)[,c(geneCol,countsCol), drop=F]
    if(all.equal(output[,1], newData[,1])){
      output = cbind(output, newData[, 2])
    } else {
      print('skipping sample: ', file)
    }
    progress_bar$tick()
  }
  
  rownames(output) = output[,1]
  output = output[,-1]
  if(!is.null(sampleNames) && length(sampleNames == ncol(output))){
    colnames(output) = sampleNames
  } else {
    colnames(output) = paste0("sample_", 1:ncol(output))
  }
  return(output)
}



# Calculates ratios between paired samples in a DESeq2ExpressionSet
# This function calculates the ratios between paired samples in a DESeq2ExpressionSet.
# 
# Arguments:
#   
    # eset: DESeq2ExpressionSet object containing the gene expression data.
    # is_log2: Logical value indicating whether the data is already in log2 scale.
    # group_names: Character vector specifying the group names.
    # groups_col: Column name in the eset object containing the group information.
    # conditions_col: Column name in the eset object containing the condition information.
    # primary_condition: Name of the primary condition for ratio calculation.
    # ref_condition: Name of the reference condition for ratio calculation.
# Returns:
    # data frame with gene names as row names and ratios between paired samples as columns.


calculate_ratios = function(eset, is_log2, group_names, groups_col, conditions_col, primary_condition, ref_condition) {
  require(progress)
  
  # Create a progress bar
  pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                        total = ncol(eset),
                        complete = "=",   # Completion bar character
                        incomplete = "-", # Incomplete bar character
                        current = ">",    # Current bar character
                        clear = FALSE,    # If TRUE, clears the bar when finish
                        width = 100)      # Width of the progress bar
  
  # Create a data frame to store the differences
  difference = data.frame(genes = rownames(eset))
  
  # Check if the data needs to be transformed to log2 scale
  if(!is_log2) {
    assay(eset) = log2(assay(eset) + 1)
  }
  
  # Iterate over each group
  for(group in group_names) {
    sub_eset = eset[, eset[[groups_col]] == group]
    
    # Check if the subset has two columns and contains the primary and reference conditions
    if(dim(sub_eset)[2] == 2 & all(sub_eset[[conditions_col]] %in% c(primary_condition, ref_condition))) {
      # Calculate the difference between primary and reference conditions
      difference[group] = assay(sub_eset[, sub_eset[[conditions_col]] == primary_condition]) - assay(sub_eset[, sub_eset[[conditions_col]] == ref_condition])
    }
    
    pb$tick()  # Update the progress bar
  }
  
  row.names(difference) = difference$genes
  difference = difference[, -1]
  
  # Transform the differences back to original scale if needed
  if(!is_log2) {
    difference = 2^difference
  }
  
  return(difference)
}



