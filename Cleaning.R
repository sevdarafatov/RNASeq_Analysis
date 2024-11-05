setwd("/mnt/DATA_4TB/projects/Anterior_Pituatary_FAUQUIER/02_analyses")

# Load the necessary library
library(dplyr)

# Load the CSV file
dt <- read.csv("DESeq2_LRT_results_without_day48.csv")

# Display the first few rows and the number of rows
head(dt)
nrow(dt)

# Remove rows with NA values
cleaned_dt <- na.omit(dt)

# Display the first few rows and the number of rows after cleaning
head(cleaned_dt)
nrow(cleaned_dt)
colnames(cleaned_dt)

# Modify the cleaned_data to extract EnsemblID and GeneName correctly
cleaned_dt <- cleaned_dt %>%
  mutate(
    EnsemblID = sub("_.*|\\|.*", "", X),  # Extract everything before the first _ or |
    
    # Conditional logic: if an underscore is present, take everything after it,
    # otherwise, take everything after the pipe
    GeneName = ifelse(grepl("_", X),
                      sub("^[^_]*_", "", X),   # Keep everything after the first underscore
                      sub("^[^|]*\\|", "", X)  # If no underscore, keep everything after the pipe
    )
  )

# Rearrange columns to put EnsemblID and GeneName first
cleaned_dt <- cleaned_dt %>%
  select(EnsemblID, GeneName, everything(), -X)  # Select EnsemblID and GeneName first, then all other columns, excluding X

write.csv(cleaned_dt, "DESeq2_LRT_results_without_day48_cleaned.csv", row.names = FALSE) 

# Filter the data based on padj threshold (e.g., padj < 0.01)
threshold <- 0.01
filtered_data <- cleaned_data %>%
  filter(padj < threshold)

# Display the first few rows of the filtered data
head(filtered_data)

# Optionally, save the filtered data back to a CSV file
write.csv(filtered_data, "DESeq2_LRT_results_without_day48_cleaned_filtered_0.01.csv", row.names = FALSE)

# Filter the data based on padj threshold (e.g., padj < 0.05)
threshold <- 0.05
filtered_data <- cleaned_data %>%
  filter(padj < threshold)

# Display the first few rows of the filtered data
head(filtered_data)

# Optionally, save the filtered data back to a CSV file
write.csv(filtered_data, "DESeq2_LRT_results_without_day48_cleaned_filtered_0.05.csv", row.names = FALSE)

