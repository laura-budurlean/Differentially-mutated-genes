# Install if needed, and then load necessary packages in R.
library(dplyr)
library(data.table)
library(writexl)

#Specify your population sizes
Population_1 <- 100
Population_2 <- 100

# Load mutation table for population 1. The format is a tab delimited text file with column headers: Gene	Variant_Type	Sample_ID
mutations1 <- fread("population1_mutations_for_Fishers_tests.txt")

# Load mutation table for population 2. The format is a tab delimited text file with column headers: Gene	Variant_Type	Sample_ID
mutations2 <- fread("population2_mutations_for_Fishers_tests.txt")

# Create vector of genes to analyze. Format is a tab delimited text file with a list of gene names to test for, and no header. 
gene_file <- "COSMIC_cancer_genes.txt" # If you are using a curated list of genes, such as known cancer drivers in a specific cancer type, you do not use adjusted P-values. If your list of genes is any gene possible, then you must use a correction. (Ie FDR)
gene_list <- read.delim(gene_file, header = FALSE, stringsAsFactors = FALSE)
genes <- gene_list$V1

#################### Start processing input files #######################

# Create empty list to store contingency tables
tables <- vector("list", length(genes))


# Loop through genes and create contingency table for each gene
for (i in seq_along(genes)) {
     gene <- genes[i]
     
     # Subset mutation tables for the current gene
     gene_mutations1 <- mutations1 %>%
          filter(Gene == gene)
     
     gene_mutations2 <- mutations2 %>%
          filter(Gene == gene)
     
  
# Calculate counts of mutated and non-mutated individuals for each population
     mutated_count1 <- sum(gene_mutations1$Variant_Type %in% c("DELETION", "INSERTION", "DUP/GAIN", "INVERSION")) # replace these with the names in the Variant_types column in your input gene list, mutations1.
     not_mutated_count1 <- Population_1 - mutated_count1
     
     mutated_count2 <- sum(gene_mutations2$Variant_Type %in% c("DELETION", "INSERTION", "DUP/GAIN", "INVERSION")) # replace these with the names in the Variant_types column in your input gene list, mutations2.
     not_mutated_count2 <- Population_2 - mutated_count2
     
 # Create contingency table
     table <- matrix(
          c(
               mutated_count1,
               not_mutated_count1,
               mutated_count2,
               not_mutated_count2
          ),
          nrow = 2,
          ncol = 2,
          dimnames = list(c("Mutated", "Not mutated"), c("Population 1", "Population 2"))
     )
     
     # Add table to list
     tables[[i]] <- table
}

# Save contingency tables to file
saveRDS(tables, "contingency_tables.rds")

# Loop through the contingency tables
for (i in seq_along(tables)) {
     table <- tables[[i]]
     
# Check if any row contains negative values. If you specified your population sizes correctly then you should not have negative values. If you have negative values the Fisher's test will not work. It usually indicates a problem with the sample counts in your input files or incorrectly specifying the number of samples you have.
     if (any(table < 0)) {
          # Print the gene name and the corresponding table
          cat("Gene:", genes[i], "\n")
          print(table)
     }
}



##################################FISHERS TEST AND CORRECTION ##############################

# Load contingency tables from file
tables <- readRDS("contingency_tables.rds")

# Perform Fisher's exact test on valid contingency tables
p_values <- sapply(tables, function(table) {
     fisher.test(table, alternative = "two.sided")$p.value}) # Note: this performs a two sided test. Adjust options as needed.

# Apply multiple testing correction if needed (e.g., Benjamini-Hochberg FDR)
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# Create a data frame with gene names, p-values, and adjusted p-values
results <- data.frame(Gene = genes, P_Value = p_values, Adjusted_P_Value = adjusted_p_values)

# Sort the results based on adjusted p-values
results <- results[order(results$Adjusted_P_Value), ]

# Print the top significant genes
top_significant <- subset(results, Adjusted_P_Value < 0.05)
print(top_significant)

# Or print all results
print(results)

# Save the results to an Excel readable file
write_xlsx(results, "Fishers_results.xlsx")



##################################CHi-SQUARE AND CORRECTION ##############################
# Instead of Fishers test, you may wish to perform chi-square test on each contingency table if applicable.
results_chi_sq <- lapply(tables, chisq.test)

# Extract p-values from the results
p_values_chi_sq <- sapply(results_chi_sq, function(res) res$p.value)

# Apply Benjamini-Hochberg correction
adjusted_p_values_chi_sq <- p.adjust(p_values_chi_sq, method = "fdr")

# Print the adjusted p-values
print(adjusted_p_values_chi_sq)

# Sort the results based on adjusted p-values
results_chi_sq <- results[order(results_chi_sq$adjusted_p_values_chi_sq), ]

# Print the top significant genes
top_significant_chi_sq <- subset(results_chi_sq, adjusted_p_values_chi_sq < 0.05)
print(top_significant_chi_sq)

# Or print all results
print(results_chi_sq)

# Save the results to an Excel readable file
write_xlsx(results, "ChiSquare_results.xlsx")
