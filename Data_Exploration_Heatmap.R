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

set.seed(123)

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]

### Log-Transformation and Filtering ###

# Log-transform the data (adding 1 to avoid log(0) issues)
log_transformed_matrix <- log1p(filtered_final_df)
head(log_transformed_matrix)

# Remove rows with zero variance
variable_genes <- apply(log_transformed_matrix, 1, var) > 0
filtered_matrix <- log_transformed_matrix[variable_genes, ]

# Define grouping variables
Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", 
               "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "WT", 
               "KI", "KI", "KI", "KI", "KI", "KI")

TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", 
               "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", 
               "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")


# Create metadata dataframe
metadata <- data.frame(
  row.names = colnames(filtered_final_df),
  Condition = Condition,
  TimePoint = TimePoint
)

# Assuming you have filtered_final_df from previous steps and metadata with batch information
metadata$Batch <- ifelse(metadata$TimePoint == "d48", "Batch2", "Batch1")


# Combine Condition and TimePoint into a single grouping variable for annotations
annotation_col <- data.frame(Condition = Condition, TimePoint = TimePoint, Batch = metadata$Batch)

# Ensure annotation_col row names match filtered_matrix column names
rownames(annotation_col) <- colnames(filtered_matrix) 

### Heatmap Generation for All Data ###

# Open PDF device to save the heatmap
pdf("heatmap_all_data_raw.pdf", width = 12, height = 10)

# Generate the heatmap for all the data
pheatmap(filtered_matrix,
         scale = "row",  # Scale genes (rows)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of RNA-Seq Data Across All Time Points",
         fontsize_row = 6,  # Adjust font size for rows
         fontsize_col = 8,  # Adjust font size for columns (samples)
         show_rownames = FALSE,  # Hide row names for large datasets
         show_colnames = TRUE,  # Show sample names (columns)
         annotation_col = annotation_col  # Add Condition and TimePoint annotations at the top
)

# Close the PDF device to finalize the file
dev.off()




### Heatmap after outlier removal and filtration ###

### Log-Transformation and Filtering ###

# Convert the filtered final data frame to a matrix
mx_filtered_final_df <- as.matrix(filtered_final_df)
head(mx_filtered_final_df)

# Log-transform the data (adding 1 to avoid log(0) issues)
log_transformed_matrix <- log1p(mx_filtered_final_df)
head(log_transformed_matrix)

# Remove rows with zero variance
row_var <- rowVars(log_transformed_matrix)  # Requires matrixStats package
variable_genes <- row_var > 0
filtered_matrix <- log_transformed_matrix[variable_genes, ]


# Define grouping variables
Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")


# Create metadata dataframe
metadata <- data.frame(
  row.names = colnames(filtered_matrix),
  Condition = Condition,
  TimePoint = TimePoint
)

ncol(filtered_final_df)
length(Condition)
length(TimePoint)


# Assuming you have filtered_final_df from previous steps and metadata with batch information
metadata$Batch <- ifelse(metadata$TimePoint == "d48", "Batch2", "Batch1")

# Combine Condition and TimePoint into a single grouping variable for annotations
annotation_col <- data.frame(Condition = Condition, TimePoint = TimePoint, Batch = metadata$Batch)

# Ensure annotation_col row names match filtered_matrix column names
rownames(annotation_col) <- colnames(filtered_matrix) 

### Heatmap Generation for All Data ###

# Open PDF device to save the heatmap
pdf("heatmap_all_data_outlier_removed_filtered_new.pdf", width = 12, height = 10)

# Generate the heatmap for all the data
pheatmap(filtered_matrix,
         scale = "row",  # Scale genes (rows)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Heatmap of RNA-Seq Data Across All Time Points",
         fontsize_row = 6,  # Adjust font size for rows
         fontsize_col = 8,  # Adjust font size for columns (samples)
         show_rownames = FALSE,  # Hide row names for large datasets
         show_colnames = TRUE,  # Show sample names (columns)
         annotation_col = annotation_col  # Add Condition and TimePoint annotations at the top
)

# Close the PDF device to finalize the file
dev.off()


