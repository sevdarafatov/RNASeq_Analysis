#### Sevda Rafatov 30.08.2024 ####

### Data exploration ###

setwd("/mnt/DATA_4TB/projects/Anterior_Pituatary_FAUQUIER/02_analyses")

# Load libraries
library(readr)
library(ggplot2)
library(RColorBrewer)  # For color palettes
library(ggrepel)
library(pheatmap)      # For heatmap generation
library(tibble)
library(DESeq2)

set.seed(123)

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]

# Remove the outlier sample 2022A115
outlier_removed_df <- final_df[, !colnames(final_df) %in% '2022A115']

# Keep only rows that have a count of at least 10 for a minimum of 3 samples
filtered_final_df <- outlier_removed_df[rowSums(outlier_removed_df >= 10) >= 3, ]

# Metadata: Condition, TimePoint, Batch
Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")
Batch <- c("Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch2", 
           "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", 
           "Batch2", "Batch2", "Batch2")

# Create metadata dataframe
metadata <- data.frame(
  row.names = colnames(filtered_final_df),
  Condition = Condition,
  TimePoint = TimePoint,
  Batch = Batch
)

### Differential expression analysis for each time point ###

# Define a function to run DESeq2 for a specific time point
run_deseq_for_timepoint <- function(timepoint) {
  # Subset metadata for the current time point
  metadata_tp <- metadata[metadata$TimePoint == timepoint, ]
  
  # Subset count data for the same samples
  counts_tp <- filtered_final_df[, rownames(metadata_tp)]
  
  # Convert Condition column to factor
  metadata_tp$Condition <- factor(metadata_tp$Condition)
  
  # Create DESeq2 dataset for the current time point
  dds_tp <- DESeqDataSetFromMatrix(countData = counts_tp,
                                   colData = metadata_tp,
                                   design = ~ Condition)
  
  # Run DESeq2
  dds_tp <- DESeq(dds_tp)
  
  # Get results for KI vs WT
  res_tp <- results(dds_tp, contrast = c("Condition", "KI", "WT"))
  
  # Return the results
  return(res_tp)
}

# Run DESeq2 for each time point
res_d27 <- run_deseq_for_timepoint("d27")
res_d48 <- run_deseq_for_timepoint("d48")
res_d75 <- run_deseq_for_timepoint("d75")
res_d105 <- run_deseq_for_timepoint("d105")

# View results summary for each time point
summary(res_d27)
summary(res_d48)
summary(res_d75)
summary(res_d105)

# Save results for each time point
write.csv(as.data.frame(res_d27), "DESeq2_results_d27.csv")
write.csv(as.data.frame(res_d48), "DESeq2_results_d48.csv")
write.csv(as.data.frame(res_d75), "DESeq2_results_d75.csv")
write.csv(as.data.frame(res_d105), "DESeq2_results_d105.csv")


# Identify top 20 DEGs by adjusted p-value
get_top_genes <- function(res, top_n = 20) {
  # Remove rows with NA adjusted p-values
  res <- res[!is.na(res$padj), ]
  
  # Sort results by adjusted p-value (smallest to largest)
  res_sorted <- res[order(res$padj), ]
  
  # Select the top N genes with the lowest adjusted p-values
  top_genes <- head(res_sorted, top_n)
  
  # Return the top genes
  return(top_genes)
}

# Get top 20 DEGs for each time point
top_genes_d27 <- get_top_genes(res_d27, top_n = 20)
top_genes_d48 <- get_top_genes(res_d48, top_n = 20)
top_genes_d75 <- get_top_genes(res_d75, top_n = 20)
top_genes_d105 <- get_top_genes(res_d105, top_n = 20)

# View the top 20 genes for day 27
top_genes_d27

write.csv(as.data.frame(top_genes_d27), "Top20_DEGs_d27.csv")
write.csv(as.data.frame(top_genes_d48), "Top20_DEGs_d48.csv")
write.csv(as.data.frame(top_genes_d75), "Top20_DEGs_d75.csv")
write.csv(as.data.frame(top_genes_d105), "Top20_DEGs_d105.csv")



### Time Course Differential Expression Analysis ###

# Ensure Condition and TimePoint are factors in the metadata
metadata$Condition <- factor(metadata$Condition)
metadata$TimePoint <- factor(metadata$TimePoint)

# Create DESeq2 dataset with the full model (including the interaction term)
dds_full <- DESeqDataSetFromMatrix(countData = filtered_final_df, 
                                   colData = metadata,
                                    design = ~ Condition + TimePoint + Condition:TimePoint)

# Run DESeq2 to estimate dispersions and size factors
dds_full <- DESeq(dds_full)

# Perform likelihood ratio test (LRT) using the reduced model (excluding the Condition:TimePoint interaction)
dds_full_lrt <- nbinomLRT(dds_full, reduced = ~ Condition + TimePoint)

# Extract results from the LRT
res_lrt <- results(dds_full_lrt)

# Check the summary of the results
summary(res_lrt)

# View top differentially expressed genes from LRT (sorted by adjusted p-value)
top_degs_lrt <- head(res_lrt[order(res_lrt$padj), ], 20)

# Save the full LRT results and top 20 DEGs
write.csv(as.data.frame(res_lrt), "DESeq2_LRT_results_all_timepoints.csv")
write.csv(as.data.frame(top_degs_lrt), "Top20_DEGs_LRT.csv")




### Time Course Differential Expression Analysis Without Day48 ###

# Ensure Condition and TimePoint are factors in the metadata
metadata$Condition <- factor(metadata$Condition)
metadata$TimePoint <- factor(metadata$TimePoint)

# Filter metadata and counts for the selected time points (d27, d75, d105)
selected_timepoints <- c("d27", "d75", "d105")
filtered_metadata <- metadata[metadata$TimePoint %in% selected_timepoints, ]

# Subset the count data based on the filtered metadata
filtered_counts <- filtered_final_df[, rownames(filtered_metadata)]

# Create DESeq2 dataset with the full model (including the interaction term)
dds_full <- DESeqDataSetFromMatrix(countData = filtered_counts, 
                                    colData = filtered_metadata,
                                    design = ~ Condition + TimePoint + Condition:TimePoint)

# Run DESeq2 to estimate dispersions and size factors
dds_full <- DESeq(dds_full)

# Perform likelihood ratio test (LRT) using the reduced model (excluding the Condition:TimePoint interaction)
dds_full_lrt <- nbinomLRT(dds_full, reduced = ~ Condition + TimePoint)

# Extract results from the LRT
res_lrt <- results(dds_full_lrt)

# Check the summary of the results
summary(res_lrt)

# View top differentially expressed genes from LRT (sorted by adjusted p-value)
top_degs_lrt <- head(res_lrt[order(res_lrt$padj), ], 10)

# Save the full LRT results and top 20 DEGs
write.csv(as.data.frame(res_lrt), "DESeq2_LRT_results_selected_timepoints.csv")
write.csv(as.data.frame(top_degs_lrt), "Top10_DEGs_LRT_selected_timepoints.csv")

