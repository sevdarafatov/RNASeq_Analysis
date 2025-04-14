
setwd("/mnt/DATA_4TB/projects/Anterior_Pituatary_FAUQUIER/02_analyses")
set.seed(123)
# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query transcripts for TBX19
tbx19_transcripts <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
                           filters = "external_gene_name",
                           values = "TBX19",
                           mart = ensembl)

print(tbx19_transcripts)


# 1. Define file paths (adjust as needed)
file1 <- "A.csv"  # First file location
file2 <- "B.csv"  # Second file location

# Read the files as CSV (since they are comma-separated)
data1 <- read_csv(file1)
data2 <- read_csv(file2)

# Preview column names to confirm "transcript_id" is present
colnames(data1)
colnames(data2)

# Merge by 'transcript_id'
merged_data <- full_join(data1, data2, by = "transcript_id")
head(merged_data)
# Save the result
write_csv(merged_data, "merged_transcript_counts.csv")  # CSV since we're using commas

# 1. Read the transcript count matrix
transcript_counts <- read_csv("merged_transcript_counts.csv")

# 2. Identify TBX19 transcripts (replace with actual Ensembl IDs if known)
# Example: Filter rows where transcript_id contains "ENST00000338863" (from your data)
tbx19_transcripts <- transcript_counts %>%
  filter(str_detect(transcript_id, "ENST00000367821|ENST00000431969|ENST00000441464|ENST00000465440")) 

tbx19_long <- tbx19_transcripts %>%
  pivot_longer(cols = starts_with("20"), names_to = "sample", values_to = "counts")

plot1 <- ggplot(tbx19_long, aes(x = sample, y = counts, fill = transcript_id)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "TBX19 Transcript Expression Across Samples",
       x = "Sample ID",
       y = "Read Counts") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_discrete(name = "Transcript ID")

ggsave("Barplot of transcripts for TBX19.pdf", plot = plot1,  width = 10, height = 8)
