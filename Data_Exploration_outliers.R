#### Sevda Rafatov 30.08.2024 ####

### Data exploration ###

setwd("/mnt/DATA_4TB/projects/Anterior_Pituatary_FAUQUIER/02_analyses")

# Install necessary packages (uncomment if not already installed)
#install.packages("ggplot2")      # For plotting
# install.packages("pheatmap")     # For visualization
#install.packages("DESeq2")       # For normalization
#install.packages("readr")        # For reading data


# Load libraries
library(readr)
library(ggplot2)
library(RColorBrewer)  # For color palettes
library(ggrepel)
library(pheatmap)      # For heatmap generation
library(tibble)

set.seed(123)

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Write the csv file #
#write_csv(merged_df, "merged_df.csv")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]
head(final_df)
nrow(final_df)


# Remove the outlier sample 2022A115 #
outlier_removed_df <- final_df[, !colnames(final_df) %in% '2022A115']
head(outlier_removed_df)
nrow(outlier_removed_df)

# Remove genes where the sum of counts across all samples is lower than 10 #
#filtered_final_df <- outlier_removed_df[rowSums(outlier_removed_df) >= 10, ]
#head(filtered_final_df)


# Keep only rows that have a count of at least 10 for a minumum of 3 samples #
filtered_final_df <- outlier_removed_df[rowSums(outlier_removed_df >= 10) >= 3, ]

# To check if the filtering worked
nrow(outlier_removed_df)   # Number of genes before filtering
nrow(filtered_final_df)    # Number of genes after filtering
min(rowSums(filtered_final_df >= 10))  # Ensure each gene has at least 3 samples with counts >= 10
head(filtered_final_df)


### Performing PCA ###

# Log-transform the data (adding 1 to avoid log(0) issues)
log_transformed_matrix <- log1p(filtered_final_df)
nrow(log_transformed_matrix)

# Remove rows with zero variance
variable_genes <- apply(log_transformed_matrix, 1, var) > 0
filtered_matrix <- log_transformed_matrix[variable_genes, ]
nrow(filtered_matrix)

#ncol(filtered_matrix)

# Perform PCA
pca <- prcomp(t(filtered_matrix), scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

write_csv(filtered_final_df, "filtered_final_df.csv")

Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")

# Combine Condition and Time Point into a single grouping variable
Group <- paste(Condition, TimePoint, sep = "_")

# Assuming you already performed PCA and have the pca_df
pca_df <- data.frame(PC1 = pca$x[,1], 
                     PC2 = pca$x[,2], 
                     Sample = rownames(pca$x),
                     Group = Group)  # Add the Group column


# Plot the PCA

# Plot the PCA result #
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Group)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +
  ggtitle("PCA of RNA-Seq Data Outlier Removed & Filtered_min_3samples_10counts") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_color_brewer(palette = "Set1")


# Save the PCA plot
ggsave("pca_plot_outlier_removed_filtered_min_3samples_10counts_06.11.2024.pdf", plot = pca_plot, width = 10, height = 8)


### Plotting PCA for batches before normalization&correction ###

# Perform PCA
pca <- prcomp(t(log_transformed_matrix), scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100



Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")
Batch <- c("Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch2", 
           "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", 
           "Batch2", "Batch2", "Batch2")

# Combine Condition and Time Point into a single grouping variable
Group <- paste(TimePoint, Condition, sep = "_")

# Assuming you already performed PCA and have the pca_df
pca_df <- data.frame(PC1 = pca$x[,1], 
                     PC2 = pca$x[,2], 
                     Sample = rownames(pca$x),
                     Group = Group)  # Add the Group column


# Plot the PCA

# Plot the PCA result #
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Batch)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +
  ggtitle("PCA of RNA-Seq Data Outlier Removed&Filtered_min_3samples_10counts") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_color_brewer(palette = "Set1")

# Save the PCA plot
ggsave("pca_plot_outlier_removed_filtered_colors_corresponding to_batches_06.11.2024.pdf", plot = pca_plot, width = 10, height = 8)




# Plot the Heatmap #

# Combine Condition and TimePoint into a single grouping variable for annotations
annotation_col <- data.frame(Condition = Condition, TimePoint = TimePoint)
# Ensure annotation_col row names match filtered_matrix column names
rownames(annotation_col) <- colnames(filtered_matrix) 

### Heatmap Generation for All Data ###

# Open PDF device to save the heatmap
pdf("heatmap_all_data_outlier_removed&filtered_min_3samples_10counts.pdf", width = 12, height = 10)

# Generate the heatmap for all the data
pheatmap(filtered_matrix,
         scale = "row",  # Scale genes (rows)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of RNA-Seq Data Outlier Removed & Filtered min_3samples_10counts",
         fontsize_row = 6,  # Adjust font size for rows
         fontsize_col = 8,  # Adjust font size for columns (samples)
         show_rownames = FALSE,  # Hide row names for large datasets
         show_colnames = TRUE,  # Show sample names (columns)
         annotation_col = annotation_col  # Add Condition and TimePoint annotations at the top
)

# Close the PDF device to finalize the file
dev.off()



# Normalization of the data #

# Install necessary packages (if not already installed)
#install.packages("DESeq2")  # For normalization
#install.packages("sva")     # For batch correction

#Load libraries

library(DESeq2)
library(sva) # For batch correction

Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")
Batch <- c("Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch2", 
           "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", 
           "Batch2", "Batch2", "Batch2")
# Combine Condition and Time Point into a single grouping variable
Group <- paste(Condition, TimePoint, sep = "_")

# Prepare a sample metadata table
# Example metadata, adjust according to your design
metadata <- data.frame(
  row.names = colnames(filtered_final_df),
  Condition = Condition,   # Your condition data
  TimePoint = TimePoint    # Your time points data
)

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered_final_df, 
                              colData = metadata, 
                              design = ~ Condition + TimePoint)


# Run the DESeq pipeline (this also normalizes the data)
dds <- DESeq(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Check the normalized counts
head(normalized_counts)



# PCA for normalized counts #
# Log-transform the corrected counts (adding 1 to avoid log(0))
log_corrected_counts2 <- log1p(normalized_counts)
head(log_corrected_counts2)

# Perform PCA on the normalized data
pca_corrected <- prcomp(t(log_corrected_counts2), scale. = TRUE)

# Create a PCA data frame
pca_df_corrected <- data.frame(PC1 = pca_corrected$x[,1], 
                               PC2 = pca_corrected$x[,2], 
                               Sample = rownames(pca_corrected$x),
                               Group = Group)

# Plot the PCA after batch correction
pca_plot3 <- ggplot(pca_df_corrected, aes(x = PC1, y = PC2, label = Sample, color = Batch)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
    xlab(paste0("PC1: ", round((pca_corrected$sdev^2 / sum(pca_corrected$sdev^2))[1] * 100, 2), "% variance")) +
    ylab(paste0("PC2: ", round((pca_corrected$sdev^2 / sum(pca_corrected$sdev^2))[2] * 100, 2), "% variance")) +
    ggtitle("PCA After Normalizaton") +
    theme_minimal()

ggsave("PCA_after_normalization_by_deseq2_colored_by_batch.pdf", plot = pca_plot3, width = 10, height = 8)


### Using SVA for batch correction ###

# Add a batch column to metadata
# Example batch column (ensure it corresponds to your actual batch data)
metadata$Batch <- c(rep("Batch1", 23), rep("Batch2", 11))  # Adjust batch labels

# Apply ComBat to correct for batch effects
corrected_counts <- ComBat_seq(as.matrix(normalized_counts), batch = metadata$Batch)

# Check the batch-corrected counts
head(corrected_counts)

# Save normalized counts
write.csv(normalized_counts, "normalized_counts.csv")

# Save batch-corrected counts
write.csv(corrected_counts, "batch_corrected_counts.csv")


# PCA After Batch Correction #

# Log-transform the corrected counts (adding 1 to avoid log(0))
log_corrected_counts <- log1p(corrected_counts)

# Perform PCA on the batch-corrected data
pca_corrected <- prcomp(t(log_corrected_counts), scale. = TRUE)

# Create a PCA data frame
pca_df_corrected <- data.frame(PC1 = pca_corrected$x[,1], 
                               PC2 = pca_corrected$x[,2], 
                               Sample = rownames(pca_corrected$x),
                               Group = Group)

# Plot the PCA after batch correction
pca_plot3 <- ggplot(pca_df_corrected, aes(x = PC1, y = PC2, label = Sample, color = Batch)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
    xlab(paste0("PC1: ", round((pca_corrected$sdev^2 / sum(pca_corrected$sdev^2))[1] * 100, 2), "% variance")) +
    ylab(paste0("PC2: ", round((pca_corrected$sdev^2 / sum(pca_corrected$sdev^2))[2] * 100, 2), "% variance")) +
    ggtitle("PCA After Normalizaton and Batch Correction") +
    theme_minimal()

ggsave("PCA_after_normalization_and_batch_correction.pdf", plot = pca_plot3, width = 10, height = 8)

ggsave("PCA_after_normalization_and_batch_correction_colors_corresponding_to_batches.pdf", plot = pca_plot3, width = 10, height = 8)


### PCA only according to the Condition ###

# Assuming you've already performed SVA and have the corrected_matrix after batch correction
# Perform PCA on the batch-corrected matrix
pca <- prcomp(t(log_corrected_counts), scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# Create a PCA dataframe, only using Condition as a grouping variable
pca_df <- data.frame(PC1 = pca$x[, 1], 
                     PC2 = pca$x[, 2], 
                     Sample = rownames(pca$x),
                     Condition = metadata$Condition)  # Add Condition column (WT/KI)

# Plot the PCA result, coloring by Condition (WT/KI)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(size = 3) +  # Plot points
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +  # Add labels
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +  # Label x-axis
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +  # Label y-axis
  ggtitle("PCA of RNA-Seq Data (Condition: WT vs KI)") +  # Title
  theme_minimal() +  # Minimal theme
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("WT" = "blue", "KI" = "red"))  # Color by Condition (WT=blue, KI=red)

# Save the PCA plot
ggsave("pca_plot_sva_batch_corrected_by_condition.pdf", plot = pca_plot, width = 10, height = 8)


### PCA only according to the Time Points ###

# Assuming you've already performed SVA and have the corrected_matrix after batch correction
# Perform PCA on the batch-corrected matrix
pca <- prcomp(t(log_corrected_counts), scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# Create a PCA dataframe, using TimePoint as a grouping variable
pca_df <- data.frame(PC1 = pca$x[, 1], 
                     PC2 = pca$x[, 2], 
                     Sample = rownames(pca$x),
                     TimePoint = metadata$TimePoint)  # Add TimePoint column

# Plot the PCA result, coloring by TimePoint
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = TimePoint)) +
  geom_point(size = 3) +  # Plot points
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +  # Add labels
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +  # Label x-axis
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +  # Label y-axis
  ggtitle("PCA of RNA-Seq Data (Colored by TimePoints)") +  # Title
  theme_minimal() +  # Minimal theme
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_color_brewer(palette = "Set1")  # Use a color palette for TimePoints

# Save the PCA plot
ggsave("pca_plot_sva_batch_corrected_by_timepoint.pdf", plot = pca_plot, width = 10, height = 8)


### Heatmap after normalization & batch correction ###

# Combine Condition and TimePoint into a single grouping variable for annotations #

annotation_col <- data.frame(Condition = Condition, TimePoint = TimePoint)

# Ensure annotation_col row names match corrected counts column names
rownames(annotation_col) <- colnames(corrected_counts)

# Proceed with Heatmap
pdf("heatmap_batch_corrected_sva.pdf", width = 12, height = 10)

pheatmap(corrected_counts,
         scale = "row",  
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of RNA-Seq Data with Batch Correction",
         fontsize_row = 6,
         fontsize_col = 8,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = annotation_col  # Add Condition and TimePoint annotations
)

dev.off()


# Using limma for batch correction #

library(limma)

# Assuming you have filtered_final_df from previous steps and metadata with batch information
metadata$Batch <- ifelse(metadata$TimePoint == "d48", "Batch2", "Batch1")

# Log-transform the data (if not done already)
log_transformed_matrix <- log1p(filtered_final_df)

# Create a design matrix for the Condition (you can also include TimePoint if needed)
design <- model.matrix(~Condition + TimePoint, data = metadata)
head(design)

# Remove the batch effect using the limma removeBatchEffect function
batch_corrected_matrix <- removeBatchEffect(log_transformed_matrix, batch = metadata$Batch)

# Now you can perform PCA on the batch-corrected data
pca_limma <- prcomp(t(batch_corrected_matrix), scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca_limma$sdev^2 / sum(pca_limma$sdev^2) * 100

# Plot PCA as previously done
pca_df <- data.frame(PC1 = pca_limma$x[, 1], 
                     PC2 = pca_limma$x[, 2], 
                     Sample = rownames(pca_limma$x),
                     Group = Group)  # Group is your Condition_TimePoint combination


pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +
  ggtitle("PCA of RNA-Seq Data with Batch Correction Limma (Mutant vs Control)") +
  theme_minimal()

# Save the PCA plot
ggsave("pca_plot_batch_corrected_by_condition_29.10.2024_.pdf", plot = pca_plot, width = 10, height = 8)




# Proceed with Heatmap

# Combine Condition and TimePoint into a single grouping variable for annotations #

annotation_col <- data.frame(Condition = Condition, TimePoint = TimePoint, Batch = Batch)

# Ensure annotation_col row names match corrected counts column names
rownames(annotation_col) <- colnames(batch_corrected_matrix)

pdf("heatmap_batch_corrected_limma_29.10.2024.pdf", width = 12, height = 10)

pheatmap(batch_corrected_matrix,
         scale = "row",  
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of RNA-Seq Data with Batch Correction Limma",
         fontsize_row = 6,
         fontsize_col = 8,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = annotation_col  # Add Condition and TimePoint annotations
)

dev.off()


### Combat_seq should be used on raw data ###

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Write the csv file #
#write_csv(merged_df, "merged_df.csv")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]
head(final_df)
nrow(final_df)


# Remove the outlier sample 2022A115 #
outlier_removed_df <- final_df[, !colnames(final_df) %in% '2022A115']
head(outlier_removed_df)
nrow(outlier_removed_df)

# Remove genes where the sum of counts across all samples is lower than 10 #
#filtered_final_df <- outlier_removed_df[rowSums(outlier_removed_df) >= 10, ]
#head(filtered_final_df)


# Keep only rows that have a count of at least 10 for a minumum of 3 samples #
filtered_final_df <- outlier_removed_df[rowSums(outlier_removed_df >= 10) >= 3, ]

# To check if the filtering worked
nrow(outlier_removed_df)   # Number of genes before filtering
nrow(filtered_final_df)    # Number of genes after filtering
min(rowSums(filtered_final_df >= 10))  # Ensure each gene has at least 3 samples with counts >= 10
head(filtered_final_df)


Batch <- c("Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch2", 
           "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", 
           "Batch2", "Batch2", "Batch2")

combat_seq_corrected <- ComBat_seq(as.matrix(filtered_final_df), batch=Batch)
head(combat_seq_corrected)

# PCA After Batch Correction #

# Log-transform the corrected counts (adding 1 to avoid log(0))
log_corrected_counts <- log1p(combat_seq_corrected)

# Perform PCA on the batch-corrected data
pca_corrected <- prcomp(t(log_corrected_counts), scale. = TRUE)

# Create a PCA data frame
pca_df_corrected <- data.frame(PC1 = pca_corrected$x[,1], 
                               PC2 = pca_corrected$x[,2], 
                               Sample = rownames(pca_corrected$x),
                               Group = Group)

# Plot the PCA after batch correction
pca_plot3 <- ggplot(pca_df_corrected, aes(x = PC1, y = PC2, label = Sample, color = Group)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
    xlab(paste0("PC1: ", round((pca_corrected$sdev^2 / sum(pca_corrected$sdev^2))[1] * 100, 2), "% variance")) +
    ylab(paste0("PC2: ", round((pca_corrected$sdev^2 / sum(pca_corrected$sdev^2))[2] * 100, 2), "% variance")) +
    ggtitle("PCA After Normalizaton and Batch Correction Combat_Seq_16.09.2024") +
    theme_minimal()

ggsave("PCA_after_normalization_and_batch_correction_groups_Combat_Seq_16092024.pdf", plot = pca_plot3, width = 10, height = 8)


# Combine Condition and TimePoint into a single grouping variable for annotations #

annotation_col <- data.frame(Condition = Condition, TimePoint = TimePoint)

# Ensure annotation_col row names match corrected counts column names
rownames(annotation_col) <- colnames(combat_seq_corrected)

# Proceed with Heatmap
pdf("heatmap_batch_CombatSeq_corrected_sva.pdf", width = 12, height = 10)

pheatmap(combat_seq_corrected,
         scale = "row",  
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of RNA-Seq Data with CombatSeq Batch Correction",
         fontsize_row = 6,
         fontsize_col = 8,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = annotation_col  # Add Condition and TimePoint annotations
)

dev.off()
