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
library(DESeq2)

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

# Keep only rows that have a count of at least 10 for a minumum of 3 samples #
filtered_final_df <- outlier_removed_df[rowSums(outlier_removed_df >= 10) >= 3, ]

# To check if the filtering worked
nrow(outlier_removed_df)   # Number of genes before filtering
nrow(filtered_final_df)    # Number of genes after filtering
min(rowSums(filtered_final_df >= 10))  # Ensure each gene has at least 3 samples with counts >= 10
head(filtered_final_df); names(filtered_final_df); dim(filtered_final_df)
                                                                            
countsData_matrix <- as.matrix(filtered_final_df)

Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")
Batch <- c("Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", 
           "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch2", 
           "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", 
           "Batch2", "Batch2", "Batch2")


# Combine Condition and Time Point into a single grouping variable
Group <- paste(Condition, TimePoint, sep = "_")

## Metadata file ##

# Prepare a sample metadata table
# Example metadata, adjust according to your design
metadata <- data.frame(
  row.names = colnames(filtered_final_df),
  Condition = Condition,   # Your condition data
  TimePoint = TimePoint    # Your time points data
)

dim(metadata)
dim(filtered_final_df)

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered_final_df, 
                              colData = metadata, 
                              design = ~ Condition + TimePoint)

# Normalization with vst #
vst <- vst(dds, blind = TRUE)
vstMatrix <- assay(vst)

head(vstMatrix)
dim(vstMatrix)


# Perform PCA on the normalized data
pca_corrected <- prcomp(t(vstMatrix), scale. = TRUE)

# Get percentage of variance explained by each principal component
pca_var <- summary(pca_corrected)$importance[2,]
pca1_var <- round(pca_var[1] * 100, 2)
pca2_var <- round(pca_var[2] * 100, 2)

# Create a data frame for PCA results
pca_df <- data.frame(
  PC1 = pca_corrected$x[,1],    # First principal component
  PC2 = pca_corrected$x[,2],    # Second principal component
  Sample = rownames(pca_corrected$x), # Sample names
  Group = Group                 # Group labels (Condition_TimePoint)
)

# Visualize PCA with ggplot2
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
  xlab(paste0("PC1: ", pca1_var, "% variance")) +
  ylab(paste0("PC2: ", pca2_var, "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()

# Save PCA plot as PDF
ggsave("PCA_plot_after_normalization_vst_according_condition_24.10.2024.pdf", plot = pca_plot, width = 10, height = 8)



# Correction of batch effect with limma #
library(limma)

# Assuming you have filtered_final_df from previous steps and metadata with batch information
metadata$Batch <- ifelse(metadata$TimePoint == "d48", "Batch2", "Batch1")


# Create a design matrix for the Condition (you can also include TimePoint if needed)
design <- model.matrix(~Condition + TimePoint, data = metadata)
head(design)

# Remove the batch effect using the limma removeBatchEffect function
batch_corrected_matrix <- removeBatchEffect(vstMatrix, batch = metadata$Batch)

# Now you can perform PCA on the batch-corrected data
pca_limma <- prcomp(t(batch_corrected_matrix), scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca_limma$sdev^2 / sum(pca_limma$sdev^2) * 100

# Plot PCA as previously done
pca_df2 <- data.frame(PC1 = pca_limma$x[, 1], 
                     PC2 = pca_limma$x[, 2], 
                     Sample = rownames(pca_limma$x),
                     Group = Group)  # Group is your Condition_TimePoint combination


# Plot PCA as previously done, using metadata$Batch for color
pca_plot <- ggplot(pca_df2, aes(x = PC1, y = PC2, label = Sample, color = metadata$TimePoint)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +
  ggtitle("PCA of RNA-Seq Data with Batch Correction Limma (Mutant vs Control)") +
  theme_minimal()

# Save the PCA plot
ggsave("pca_plot_batch_corrected_limma_time_points_09.10.2024.pdf", plot = pca_plot, width = 10, height = 8)
