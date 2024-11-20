# Load necessary libraries
library(ggplot2)

# Load the CSV file
df_genes <- read.csv("extracted_tpm_genes_of_interest.csv")
colnames(df_genes)[-1] <- gsub("^X", "", colnames(df_genes)[-1])
head(df_genes)

# Define the samples, time points, and conditions
samples <- c("2023A39", "2023A42", "2023A35", "2023A36", "2023A48", "2023A41", "2023A49", "2023A46",
             "2023A37", "2023A56", "2023A38", "2023A47", "2023A45", "2023A43", "2023A53", "2023A40",
             "2023A31", "2023A51", "2023A32", "2023A52", "2023A34", "2023A50", "2023A44", "2022A112",
             "2022A114", "2022A118", "2022A119", "2022A117", "2022A113", "2022A120", "2022A116",
             "2022A111", "2022A110", "2022A121")

Conditions <- c("WT", "WT", "KI", "KI", "WT", "WT", "WT", "KI", "KI", "KI", "KI",
                "WT", "KI", "KI", "KI", "WT", "WT", "KI", "WT", "KI", "WT", "WT",
                "KI", "WT", "WT", "KI", "KI", "KI", "WT", "KI", "KI", "WT", "WT", "KI")

Time_Point <- c("d75", "d75", "d27", "d27", "d105", "d75", "d105", "d75", "d27", "d105", "d27",
                "d105", "d75", "d75", "d105", "d75", "d27", "d105", "d27", "d105", "d27",
                "d105", "d75", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48", "d48",
                "d48", "d48")

# Open a PDF device to save all plots in one file
pdf("Boxplots_of_TPM_Values_by_Time_Point_and_Condition.pdf", width = 10, height = 8)

# Iterate over each gene to create the plots
for (i in 1:nrow(df_genes)) {
  gene_name <- df_genes$GeneName[i]
  values <- as.numeric(df_genes[i, -c(1, 2)])  # Extract TPM values (excluding the first two columns)
  
  # Create a data frame for plotting
  df <- data.frame(Sample = samples, TimePoint = Time_Point, Condition = Conditions, TPM = values)
  
  # Ensure that TimePoint is a factor with the correct order
  df$TimePoint <- factor(df$TimePoint, levels = c("d27", "d48", "d75", "d105"))
  
  # Ensure that Condition is a factor with the correct order (WT first)
  df$Condition <- factor(df$Condition, levels = c("WT", "KI"))
  
  # Create a new factor for the x-axis in the desired order
  df$Group <- factor(paste(df$TimePoint, df$Condition, sep = " "),
                     levels = c("d27 WT", "d27 KI", "d48 WT", "d48 KI", 
                                "d75 WT", "d75 KI", "d105 WT", "d105 KI"))
  
  # Plot the data with grouping by the new Group variable
  plot <- ggplot(df, aes(x = Group, y = TPM, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +  # Hide outliers in the boxplot
    geom_jitter(aes(color = Condition), position = position_jitter(width = 0.05), alpha = 0.7) +  # Add jittered points
    labs(title = paste("Boxplot of TPM Values for", gene_name, "by Time Point and Condition"), 
         x = "Time Point and Condition", 
         y = "TPM") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = c("WT" = "blue", "KI" = "red")) +  # Customize point colors
    scale_fill_manual(values = c("WT" = "#00aeff", "KI" = "#f08b8b"))  # Customize fill colors
  
  # Print the plot to the PDF
  print(plot)
}

# Close the PDF device
dev.off()
