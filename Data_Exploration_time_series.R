set.seed(123)

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]

# Log-transform the data
log_transformed_matrix <- log1p(final_df)

# Remove rows with zero variance
variable_genes <- apply(log_transformed_matrix, 1, var) > 0
filtered_matrix <- log_transformed_matrix[variable_genes, ]

# Define time points and conditions
Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")

# Create a data frame with the TimePoint and Condition information
pca_data <- data.frame(t(filtered_matrix), TimePoint = TimePoint, Condition = Condition)

# Filter data for day 48
day_48_data <- pca_data[pca_data$TimePoint == "d48", ]
head(day_48_data)
# Remove columns with zero variance (i.e., constant columns)
day_48_data <- day_48_data[, apply(day_48_data[, -c(ncol(day_48_data)-1, ncol(day_48_data))], 2, var) > 0]

# Perform PCA on day 48 data
pca_day_48 <- prcomp(day_48_data[, -c(ncol(day_48_data)-1, ncol(day_48_data))], scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca_day_48$sdev^2 / sum(pca_day_48$sdev^2) * 100

# Create PCA data frame
pca_df <- data.frame(PC1 = pca_day_48$x[,1], 
                     PC2 = pca_day_48$x[,2], 
                     Sample = rownames(pca_day_48$x), 
                     Condition = day_48_data$Condition)

# Plot the PCA
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(size = 3) +  # Slightly increase point size
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +  # Use ggrepel for better text positioning
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +  # X-axis label with variance explained
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +  # Y-axis label with variance explained
  ggtitle("PCA of RNA-Seq Data - Day 48 (KI vs WT)") +  # Title with day and condition information
  theme_minimal() +  # Minimal theme for clean look
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size for better readability
    axis.text = element_text(size = 12),  # Increase axis text size
    plot.title = element_text(size = 16),  # Increase plot title size
    legend.title = element_text(size = 12),  # Legend title size
    legend.text = element_text(size = 10)  # Legend text size
  ) +
  scale_color_manual(values = c("WT" = "blue", "KI" = "red")) +  # Set colors for conditions
  coord_cartesian(xlim = c(min(pca_df$PC1) - 5, max(pca_df$PC1) + 5), ylim = c(min(pca_df$PC2) - 5, max(pca_df$PC2) + 5))  # Adjust plot boundaries

# Save the plot
ggsave("pca_plot_day_48_KI_vs_WT.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave("pca_plot_day_48_KI_vs_WT.pdf", plot = pca_plot, width = 10, height = 8, units = "in")


#For day 27#
set.seed(123)

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]

# Log-transform the data
log_transformed_matrix <- log1p(final_df)

# Remove rows with zero variance
variable_genes <- apply(log_transformed_matrix, 1, var) > 0
filtered_matrix <- log_transformed_matrix[variable_genes, ]

# Define time points and conditions
Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")

# Create a data frame with the TimePoint and Condition information
pca_data <- data.frame(t(filtered_matrix), TimePoint = TimePoint, Condition = Condition)

# Filter data for day 48
day_27_data <- pca_data[pca_data$TimePoint == "d27", ]
head(day_27_data)
# Remove columns with zero variance (i.e., constant columns)
day_27_data <- day_27_data[, apply(day_27_data[, -c(ncol(day_27_data)-1, ncol(day_27_data))], 2, var) > 0]

# Perform PCA on day 27 data
pca_day_27 <- prcomp(day_27_data[, -c(ncol(day_27_data)-1, ncol(day_27_data))], scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca_day_27$sdev^2 / sum(pca_day_27$sdev^2) * 100

# Create PCA data frame
pca_df <- data.frame(PC1 = pca_day_27$x[,1], 
                     PC2 = pca_day_27$x[,2], 
                     Sample = rownames(pca_day_27$x), 
                     Condition = day_27_data$Condition)

# Plot the PCA
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(size = 3) +  # Slightly increase point size
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +  # Use ggrepel for better text positioning
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +  # X-axis label with variance explained
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +  # Y-axis label with variance explained
  ggtitle("PCA of RNA-Seq Data - Day 27 (KI vs WT)") +  # Title with day and condition information
  theme_minimal() +  # Minimal theme for clean look
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size for better readability
    axis.text = element_text(size = 12),  # Increase axis text size
    plot.title = element_text(size = 16),  # Increase plot title size
    legend.title = element_text(size = 12),  # Legend title size
    legend.text = element_text(size = 10)  # Legend text size
  ) +
  scale_color_manual(values = c("WT" = "blue", "KI" = "red")) +  # Set colors for conditions
  coord_cartesian(xlim = c(min(pca_df$PC1) - 5, max(pca_df$PC1) + 5), ylim = c(min(pca_df$PC2) - 5, max(pca_df$PC2) + 5))  # Adjust plot boundaries

# Save the plot
ggsave("pca_plot_day_27_KI_vs_WT.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave("pca_plot_day_27_KI_vs_WT.pdf", plot = pca_plot, width = 10, height = 8, units = "in")


# For day 75 #
set.seed(123)

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]

# Log-transform the data
log_transformed_matrix <- log1p(final_df)

# Remove rows with zero variance
variable_genes <- apply(log_transformed_matrix, 1, var) > 0
filtered_matrix <- log_transformed_matrix[variable_genes, ]

# Define time points and conditions
Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")

# Create a data frame with the TimePoint and Condition information
pca_data <- data.frame(t(filtered_matrix), TimePoint = TimePoint, Condition = Condition)

# Filter data for day 75
day_75_data <- pca_data[pca_data$TimePoint == "d75", ]
head(day_75_data)
# Remove columns with zero variance (i.e., constant columns)
day_75_data <- day_75_data[, apply(day_75_data[, -c(ncol(day_75_data)-1, ncol(day_75_data))], 2, var) > 0]

# Perform PCA on day 75 data
pca_day_75 <- prcomp(day_75_data[, -c(ncol(day_75_data)-1, ncol(day_75_data))], scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca_day_75$sdev^2 / sum(pca_day_75$sdev^2) * 100

# Create PCA data frame
pca_df <- data.frame(PC1 = pca_day_75$x[,1], 
                     PC2 = pca_day_75$x[,2], 
                     Sample = rownames(pca_day_75$x), 
                     Condition = day_75_data$Condition)

# Plot the PCA
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(size = 3) +  # Slightly increase point size
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +  # Use ggrepel for better text positioning
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +  # X-axis label with variance explained
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +  # Y-axis label with variance explained
  ggtitle("PCA of RNA-Seq Data - Day 75 (KI vs WT)") +  # Title with day and condition information
  theme_minimal() +  # Minimal theme for clean look
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size for better readability
    axis.text = element_text(size = 12),  # Increase axis text size
    plot.title = element_text(size = 16),  # Increase plot title size
    legend.title = element_text(size = 12),  # Legend title size
    legend.text = element_text(size = 10)  # Legend text size
  ) +
  scale_color_manual(values = c("WT" = "blue", "KI" = "red")) +  # Set colors for conditions
  coord_cartesian(xlim = c(min(pca_df$PC1) - 5, max(pca_df$PC1) + 5), ylim = c(min(pca_df$PC2) - 5, max(pca_df$PC2) + 5))  # Adjust plot boundaries

# Save the plot
ggsave("pca_plot_day_75_KI_vs_WT.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave("pca_plot_day_75_KI_vs_WT.pdf", plot = pca_plot, width = 10, height = 8, units = "in")


# For day 105 #

set.seed(123)

# Load the gene count matrices
gene_count_matrix1 <- read_csv("gene_count_matrix_d27_75_105.csv", show_col_types = FALSE)
gene_count_matrix2 <- read_csv("gene_count_matrix_d48.csv", show_col_types = FALSE)

# Merge the data frames by gene_id
merged_df <- merge(gene_count_matrix1, gene_count_matrix2, by = "gene_id")

# Set row names to GeneID and remove the GeneID column
rownames(merged_df) <- merged_df$gene_id
final_df <- merged_df[, -1]

# Log-transform the data
log_transformed_matrix <- log1p(final_df)

# Remove rows with zero variance
variable_genes <- apply(log_transformed_matrix, 1, var) > 0
filtered_matrix <- log_transformed_matrix[variable_genes, ]

# Define time points and conditions
Condition <- c("WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "WT", "WT", "WT", "WT", "WT", "WT", "KI", "KI", "KI", "KI", "KI", "KI")
TimePoint <- c("d27", "d27", "d27", "d27", "d27", "d27", "d27", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d75", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d105", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48")

# Create a data frame with the TimePoint and Condition information
pca_data <- data.frame(t(filtered_matrix), TimePoint = TimePoint, Condition = Condition)

# Filter data for day 105
day_105_data <- pca_data[pca_data$TimePoint == "d105", ]
head(day_105_data)
# Remove columns with zero variance (i.e., constant columns)
day_105_data <- day_105_data[, apply(day_105_data[, -c(ncol(day_105_data)-1, ncol(day_105_data))], 2, var) > 0]

# Perform PCA on day 105 data
pca_day_105 <- prcomp(day_105_data[, -c(ncol(day_105_data)-1, ncol(day_105_data))], scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- pca_day_105$sdev^2 / sum(pca_day_105$sdev^2) * 100

# Create PCA data frame
pca_df <- data.frame(PC1 = pca_day_105$x[,1], 
                     PC2 = pca_day_105$x[,2], 
                     Sample = rownames(pca_day_105$x), 
                     Condition = day_105_data$Condition)

# Plot the PCA
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(size = 3) +  # Slightly increase point size
  geom_text_repel(size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = 10) +  # Use ggrepel for better text positioning
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +  # X-axis label with variance explained
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +  # Y-axis label with variance explained
  ggtitle("PCA of RNA-Seq Data - Day 105 (KI vs WT)") +  # Title with day and condition information
  theme_minimal() +  # Minimal theme for clean look
  theme(
    axis.title = element_text(size = 14),  # Increase axis title size for better readability
    axis.text = element_text(size = 12),  # Increase axis text size
    plot.title = element_text(size = 16),  # Increase plot title size
    legend.title = element_text(size = 12),  # Legend title size
    legend.text = element_text(size = 10)  # Legend text size
  ) +
  scale_color_manual(values = c("WT" = "blue", "KI" = "red")) +  # Set colors for conditions
  coord_cartesian(xlim = c(min(pca_df$PC1) - 5, max(pca_df$PC1) + 5), ylim = c(min(pca_df$PC2) - 5, max(pca_df$PC2) + 5))  # Adjust plot boundaries

# Save the plot
ggsave("pca_plot_day_105_KI_vs_WT.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave("pca_plot_day_105_KI_vs_WT.pdf", plot = pca_plot, width = 10, height = 8, units = "in")
