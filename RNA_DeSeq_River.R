---
  title: "DESEQ River"
author: "Steven Carrell"
date: "2024-02-24"
output: pdf_document
---

  
setwd("F:/River")

##################################################
#########################
#########################
################################################## Conditions vs controls 15 - 2Fold
#########################
#########################
##################################################
###############Group A is control (both genders) /compare all against A
#change labels for lengths/bmi outputs for morphometrics
#run 2way anova for morphometrics
#remove sex as factor first
#secondary analysis on sex comparisons later
# look for coexposure effect in groups E/F
#cluster heatmap to confirm outliers
#median across treatments and reheatma


#Load required libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("apeglm")
#install.packages("DESeq2")
#BiocManager::install("DESeq2")
# Load necessary libraries if not already loaded
#if (!requireNamespace("ggplot2", quietly = TRUE)) {
#  install.packages("ggplot2")
#}
#if (!requireNamespace("ggfortify", quietly = TRUE)) {
#  install.packages("ggfortify")
#}
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(tximport)
library(DESeq2)
library(readr)
library(Select)
library(EnhancedVolcano)
library(apeglm)
library(dplyr)
library(DESeq2)

#############Import Data###############

# import count data
count_data <- read.table("F:/River/summary_sum2.txt", header = TRUE, row.names = 1, sep = "\t")



#import metadata 
#metadata <- read.delim("metadata.txt")
metadata2 <- read.delim("metadata2.txt")
View(metadata2)


# Subset metadata to only include conditions A, B, and C
metadata <- metadata2[metadata2$Condition %in% c("A", "B", "C"), ]

# Extract sample names from metadata
sample_names <- metadata$Samples

# Subset count data based on sample names
count_data_subset <- count_data[, sample_names]
print(dim(count_data_subset))

# Convert Desc column to factor
metadata$Desc <- factor(metadata$Desc)

# Create a new Group column to differentiate Controls from Non-Controls
metadata$Group <- ifelse(metadata$Desc == "Control", "Untreated", "Treated")
metadata$Group <- factor(metadata$Group)  # Convert Group to factor

# Print current column names of count data
print(colnames(count_data_subset))

# Print sample and description columns from metadata
print(metadata$Samples)
print(metadata$Desc)

# Create a named vector for renaming columns
#rename_vector <- setNames(metadata$Desc, metadata$Samples)

# Renaming the columns
#colnames(count_data_subset) <- rename_vector[colnames(count_data_subset)]

# Check new column names
print(colnames(count_data_subset))

# Remove the pattern '.#' from column names
#colnames(count_data_subset) <- gsub("\\.[0-9]+", "", colnames(count_data_subset))

# Print new column names to verify changes
print(colnames(count_data_subset))

# Find unique column names
unique_cols <- unique(colnames(count_data_subset))
print(unique_cols)


############################################
# # Compute median data for unique column names
# median_data <- sapply(unique_cols, function(col) {
#   # Ensure 'col' is a character string
#   col <- as.character(col)
#   
#   # Columns that match the current unique name
#   cols_to_median <- which(colnames(count_data_subset) == col)
#   
#   # Apply median to each row for these columns, ensuring matrix format even with one column
#   apply(count_data_subset[, cols_to_median, drop = FALSE], 1, median)
# })
# 
# # Convert median_data to a data frame
# median_data_df <- as.data.frame(median_data)
# 
# # Round median_data_df to nearest integer
# median_data_df <- round(median_data_df)
# 
# # Convert median_data_df to matrix
# countData <- as.matrix(median_data_df)
# 
# # Ensure column names are set correctly
# colnames(countData) <- colnames(median_data_df)
###########################################


# Ensure dplyr is loaded to use distinct
unique_metadata <- distinct(metadata, Samples, .keep_all = TRUE)
unique_metadata <- unique_metadata[, !colnames(unique_metadata) %in% c("Samples", "Condition")]

# Set the 'Desc' column as row names
rownames(metadata) <- metadata$Samples

# Exclude the 'Condition' column from metadata
metadata <- metadata[, !colnames(metadata) %in% c("Condition")]

# Rename 'Desc' to 'Condition' using base R
names(metadata)[names(metadata) == "Desc"] <- "Condition"
#metadata3 <- metadata
#metadata3 <- metadata3 %>% distinct(Condition, .keep_all = TRUE)

#rownames(metadata3) <- metadata3$Condition

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data_subset,
                              colData = metadata,
                              design = ~ Condition)

# Relevel 'Condition' to set 'Control' as the reference level
dds$Group <- relevel(dds$Group, ref = "Untreated")

# Run DESeq
dds <- DESeq(dds)

# Initialize a list to store results
results_list <- list()

levels_to_test <- setdiff(levels(metadata$Condition), "Control")
print(levels_to_test)
# Compute results for each level (sample) vs control and save to CSV
for (level in levels_to_test) {
  res <- results(dds, contrast = c("Condition", level, "Control"))
  results_list[[level]] <- res
  
  # Remove rows with NA values in padj or log2FoldChange
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  
  # Save genes, fold change, and p-value for each sample against control to CSV
  write.csv(res[, c("log2FoldChange", "padj")], paste0("DEGs_", level, "_vs_control.csv"), row.names = TRUE)
}

# Function to count DEGs at different fold changes
count_DEGs <- function(res, fold_change) {
  sum(res$padj < 0.05 & (abs(res$log2FoldChange) >= fold_change), na.rm = TRUE)
}

# Initialize a data frame to store DEG counts
deg_counts <- data.frame(Sample = character(),
                         FoldChange1 = integer(),
                         FoldChange15 = integer(),
                         FoldChange2 = integer(),
                         stringsAsFactors = FALSE)

# Loop over results, count DEGs, and save gene info to CSV
for (sample in names(results_list)) {
  res <- results_list[[sample]]
  
  # Remove rows with NA values in padj or log2FoldChange
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  
  # Count DEGs for different fold changes
  deg_counts <- rbind(deg_counts, data.frame(Sample = sample,
                                             FoldChange1 = count_DEGs(res, 1),
                                             FoldChange15 = count_DEGs(res, 1.5),
                                             FoldChange2 = count_DEGs(res, 2)))
  
  # Filter DEGs (padj < 0.05)
  deg_res <- res[res$padj < 0.05, ]
  
}

# Print the DEG counts
print(deg_counts)
write.csv(deg_counts, "deg_counts.csv", row.names = FALSE)

# Calculate percentage change between FoldChange1 and FoldChange15
deg_counts$FC1to15 <- (deg_counts$FoldChange15 / deg_counts$FoldChange1) * 100
deg_counts$FC15to2 <- (deg_counts$FoldChange2 / deg_counts$FoldChange15) * 100


# Print the updated deg_counts
print(deg_counts)

# Calculate percentage change between FoldChange1 and FoldChange15
deg_counts$FC1to15 <- (deg_counts$FoldChange15 / deg_counts$FoldChange1) * 100
deg_counts$FC15to2 <- (deg_counts$FoldChange2 / deg_counts$FoldChange15) * 100

# Initialize new table
new_table <- data.frame(Sample = character(), Value = numeric(), Source = character(), stringsAsFactors = FALSE)

# Function to determine the value and its source to add to the new table
determine_value <- function(row) {
  if (!is.na(row$FC1to15) && row$FC1to15 < 67) {
    return(data.frame(Sample = row$Sample, Value = row$FoldChange1, Source = "FoldChange1"))
  } else if (!is.na(row$FC1to15) && row$FC1to15 >= 67) {
    if (!is.na(row$FC15to2) && row$FC15to2 < 67) {
      return(data.frame(Sample = row$Sample, Value = row$FoldChange15, Source = "FoldChange15"))
    } else if (!is.na(row$FC15to2) && row$FC15to2 >= 67) {
      return(data.frame(Sample = row$Sample, Value = row$FoldChange2, Source = "FoldChange2"))
    } else {
      return(data.frame(Sample = row$Sample, Value = NA, Source = NA))
    }
  } else {
    return(data.frame(Sample = row$Sample, Value = NA, Source = NA))
  }
}

# Loop over each row in deg_counts and apply the function
for (i in 1:nrow(deg_counts)) {
  new_row <- determine_value(deg_counts[i, ])
  new_table <- rbind(new_table, new_row)
}

# Print the new table
print(new_table)
print(deg_counts)
# Merge new_table columns into deg_counts
deg_counts <- merge(deg_counts, new_table, by = "Sample")

# Print the updated deg_counts
print(deg_counts)

# Function to count DEGs with padj < 0.05
count_DEGs2 <- function(res) {
  sum(res$padj < 0.05, na.rm = TRUE)
}

# Initialize a data frame to store DEG counts
deg_counts2 <- data.frame(Sample = character(),
                          DEG_Count_padj_0.05 = integer(),
                          stringsAsFactors = FALSE)

# Loop over results and count DEGs
for (sample in names(results_list)) {
  res <- results_list[[sample]]
  
  # Count DEGs with padj < 0.05
  deg_count2 <- count_DEGs2(res)
  
  # Store the results in the data frame
  deg_counts2 <- rbind(deg_counts2, data.frame(Sample = sample,
                                               DEG_Count_padj_0.05 = deg_count2))
}

# Print the DEG counts
print(deg_counts2)

# Function to count DEGs with padj < 0.01
count_DEGs3 <- function(res) {
  sum(res$padj < 0.01, na.rm = TRUE)
}

# Initialize a data frame to store DEG counts
deg_counts3 <- data.frame(Sample = character(),
                          DEG_Count_padj_0.01 = integer(),
                          stringsAsFactors = FALSE)

# Loop over results and count DEGs
for (sample in names(results_list)) {
  res <- results_list[[sample]]
  
  # Count DEGs with padj < 0.05
  deg_count3 <- count_DEGs3(res)
  
  # Store the results in the data frame
  deg_counts3 <- rbind(deg_counts3, data.frame(Sample = sample,
                                               DEG_Count_padj_0.01 = deg_count3))
}

# Print the DEG counts
print(deg_counts2)
print(deg_counts3)

# # Secondary analysis for samples with FoldChange1 < 10
# secondary_analysis <- function(sample, res) {
#   genes_p01 <- rownames(res)[res$padj < 0.05] #Check numeric as character
#   count_genes_p01 <- length(genes_p01)
#   return(count_genes_p01)
# }
# 
# # Initialize a new column for gene count with p-value < 0.01
# deg_counts$GeneCount_p01 <- NA
# 
# 
# # Loop over each sample with FoldChange1 < 10
# for (i in 1:nrow(deg_counts)) {
#    {
#     sample <- deg_counts$Sample[i]
#     res <- results_list[[sample]]
#     gene_count <- secondary_analysis(sample, res)
#     deg_counts$GeneCount_p01[i] <- gene_count
#   }
# }

# Print the updated deg_counts with GeneCount_p01
print(deg_counts)

# Import concentrations and merge concentrations with deg_counts by Sample
concentrations <- read.csv("RIVERConcentrationdata.csv")
final_data <- merge(deg_counts, deg_counts2, by = "Sample")
final_data <- merge(final_data, deg_counts3, by = "Sample")
final_data <- merge(final_data, concentrations, by = "Sample")

final_data <- final_data %>% select(-c(X, X.1))

# Print the final merged data
print(final_data)

write.csv(final_data, "final_data.csv", row.names = FALSE)

# Ensure metadata has row names matching the sample names
rownames(metadata) <- metadata$Samples

# List to store results for each sample
results_list <- list()

# Compute results for each level (sample) vs control and save to CSV
for (level in levels_to_test) {
  res <- results(dds, contrast = c("Condition", level, "Control"))
  results_list[[level]] <- res
  
  # Remove rows with NA values in padj or log2FoldChange
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  
  # Save genes, fold change, and p-value for each sample against control to CSV
  write.csv(res[, c("log2FoldChange", "padj")], paste0("DEGs_", level, "_vs_control.csv"), row.names = TRUE)
}

# Function to read and combine CSV files for each sample
combine_csv_files <- function(levels_to_test) {
  combined_data <- data.frame()
  
  for (level in levels_to_test) {
    file_name <- paste0("DEGs_", level, "_vs_control.csv")
    if (file.exists(file_name)) {
      temp_data <- read.csv(file_name, row.names = 1)
      temp_data$Sample <- level
      combined_data <- rbind(combined_data, temp_data)
    }
  }
  
  return(combined_data)
}

# Combine all CSV files into one for all samples
combined_data <- combine_csv_files(levels_to_test)

# Save the combined data to a new CSV file
write.csv(combined_data, "Combined_DEGs_all_samples_vs_control.csv", row.names = TRUE)

# Print the combined data to verify
print(head(combined_data))

# Filter combined data by padj <= 0.05
filtered_combined_data <- combined_data[combined_data$padj <= 0.05, ]

# Save the filtered combined data to a new CSV file
#write.csv(filtered_combined_data, "Filtered_Combined_DEGs_all_samples_vs_control.csv", row.names = TRUE)

# Print the filtered combined data to verify
print(head(filtered_combined_data))

# Further filter combined data by log2FoldChange >= 1 or <= -1
filtered_combined_data <- filtered_combined_data[abs(filtered_combined_data$log2FoldChange) >= 1, ]

# Save the further filtered combined data to a new CSV file
write.csv(filtered_combined_data, "Filtered_Combined_DEGs_log2FoldChange_all_samples_vs_control.csv", row.names = TRUE)

# Print the further filtered combined data to verify
print(head(filtered_combined_data))

# Copy the rownames to a new column named Gene in filtered_combined_data
filtered_combined_data$Gene <- rownames(filtered_combined_data)

# Find duplicate row names and their counts
duplicate_counts <- as.data.frame(table(filtered_combined_data$Gene))
duplicate_counts <- duplicate_counts[duplicate_counts$Freq > 1, ]

# Rename columns for clarity
colnames(duplicate_counts) <- c("Gene", "Count")

# Save the duplicate counts data to a new CSV file
write.csv(duplicate_counts, "Duplicate_Row_Names_Count.csv", row.names = FALSE)

# Print the duplicate counts data to verify
print(head(duplicate_counts))


setwd("F:/River/DEGCSV")
#install.packages("DBI")
#install.packages("RSQLite")

library(DBI)
library(RSQLite)

ensembl_data <- read.table("ensembl_1_to_1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create a SQLite database file
con <- dbConnect(RSQLite::SQLite(), "ensembl_database.sqlite")

# Write the data to the database
dbWriteTable(con, "ensembl_1_to_1", ensembl_data, overwrite = TRUE)

# List tables in the database to verify
dbListTables(con)

# Close the connection when done
dbDisconnect(con)

csv_files <- list.files(pattern = "\\.csv$")

for (file in csv_files) {
  # Read the CSV file
  data <- read.csv(file, stringsAsFactors = FALSE)
  
  # Create a SQLite database with the same name as the CSV file (without the extension)
  db_name <- paste0(tools::file_path_sans_ext(file), ".sqlite")
  con <- dbConnect(RSQLite::SQLite(), db_name)
  
  # Match and append columns
  matched_data <- merge(data, ensembl_data[, c(1, 3, 4)], by.x = names(data)[1], by.y = names(ensembl_data)[4], all.x = TRUE)
  
  # Rename the new columns
  names(matched_data)[ncol(matched_data)-1] <- "GeneID"
  names(matched_data)[ncol(matched_data)] <- "GeneName"
  
  # Write the modified data back to the database
  table_name <- tools::file_path_sans_ext(file)
  dbWriteTable(con, table_name, matched_data, overwrite = TRUE)
  
  # Close the connection
  dbDisconnect(con)
  
  cat("Processed and updated database for:", file, "\n")
}

for (file in csv_files) {
  # Create the SQLite database name
  db_name <- paste0(tools::file_path_sans_ext(file), ".sqlite")
  
  # Connect to the SQLite database
  con <- dbConnect(RSQLite::SQLite(), db_name)
  
  # Get the table name (assuming it's the same as the CSV file name without extension)
  table_name <- tools::file_path_sans_ext(file)
  
  # Read the table from the database
  data <- dbReadTable(con, table_name)
  
  # Define the new CSV file name
  new_csv_name <- paste0(tools::file_path_sans_ext(file), "_updated.csv")
  
  # Write the data to a new CSV file
  write.csv(data, new_csv_name, row.names = FALSE)
  
  # Close the connection
  dbDisconnect(con)
  
  cat("Exported:", new_csv_name, "\n")
}

setwd("F:/River/DEGCSV")

library(gprofiler2)
library(DBI)
library(RSQLite)
library(tidyr)

ensembl_data <- read.table("ensembl_1_to_1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create a SQLite database file
con <- dbConnect(RSQLite::SQLite(), "ensembl_database.sqlite")

# Write the data to the database
dbWriteTable(con, "ensembl_1_to_1", ensembl_data, overwrite = TRUE)

# List tables in the database to verify
dbListTables(con)

# Close the connection when done
dbDisconnect(con)

csv_files <- list.files(pattern = "_vs_control\\.csv$")

for (file in csv_files) {
  # Read the CSV file
  data <- read.csv(file, stringsAsFactors = FALSE)
  
  # Create a SQLite database with the same name as the CSV file (without the extension)
  db_name <- paste0(tools::file_path_sans_ext(file), ".sqlite")
  con <- dbConnect(RSQLite::SQLite(), db_name)
  
  # Match and append columns
  matched_data <- merge(data, ensembl_data[, c(1, 3, 4)], by.x = names(data)[1], by.y = names(ensembl_data)[4], all.x = TRUE)
  
  # Rename the new columns
  names(matched_data)[ncol(matched_data)-1] <- "GeneID"
  names(matched_data)[ncol(matched_data)] <- "GeneName"
  
  # Write the modified data back to the database
  table_name <- tools::file_path_sans_ext(file)
  dbWriteTable(con, table_name, matched_data, overwrite = TRUE)
  
  # Close the connection
  dbDisconnect(con)
  
  cat("Processed and updated database for:", file, "\n")
}

for (file in csv_files) {
  # Create the SQLite database name
  db_name <- paste0(tools::file_path_sans_ext(file), ".sqlite")
  
  # Connect to the SQLite database
  con <- dbConnect(RSQLite::SQLite(), db_name)
  
  # Get the table name (assuming it's the same as the CSV file name without extension)
  table_name <- tools::file_path_sans_ext(file)
  
  # Read the table from the database
  data <- dbReadTable(con, table_name)
  
  # Define the new CSV file name
  new_csv_name <- paste0(tools::file_path_sans_ext(file), "_updated.csv")
  
  # Write the data to a new CSV file
  write.csv(data, new_csv_name, row.names = FALSE)
  
  # Close the connection
  dbDisconnect(con)
  
  cat("Exported:", new_csv_name, "\n")
}

# Run gprofiler on each updated CSV file using Zebrafish (Danio rerio) organism
updated_csv_files <- list.files(pattern = "_vs_control_updated\\.csv$")

for (updated_file in updated_csv_files) {
  # Read the updated CSV file
  updated_data <- read.csv(updated_file, stringsAsFactors = FALSE)
  
  # Filter based on padj <= 0.05 and |log2FoldChange| > 1
  filtered_data <- subset(updated_data, padj <= 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
  
  # Extract the genes from the first column of the filtered data
  genes <- filtered_data[[1]]  # First column
  
  # Run gprofiler for gene enrichment analysis using Zebrafish (Danio rerio)
  gp_results <- gost(query = genes, organism = "drerio")
  
  # Check if results exist
  if (!is.null(gp_results$result)) {
    # Flatten any list columns
    flat_results <- gp_results$result %>%
      unnest(cols = where(is.list), keep_empty = TRUE)
    
    # Export the flattened results to a CSV file
    gprofiler_csv <- paste0(tools::file_path_sans_ext(updated_file), "_gprofiler_results.csv")
    write.csv(flat_results, gprofiler_csv, row.names = FALSE)
    cat("gProfiler results saved for:", updated_file, "\n")
  } else {
    cat("No gProfiler results for:", updated_file, "\n")
  }
}

setwd("F:/River")
# Function to get the top 30 DEGs
get_top_30_DEGs <- function(res) {
  res_ordered <- res[order(res$padj, na.last = NA), ]
  top_30 <- head(res_ordered, 30)
  return(top_30)
}

for (sample in names(results_list)) {
  res <- results_list[[sample]]
  top_30 <- get_top_30_DEGs(res)
  
  # Create a bar plot for the top 30 DEGs
  plot <- ggplot(top_30, aes(x = reorder(rownames(top_30), log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top 30 DEGs in", sample),
         x = "Genes",
         y = "Log2 Fold Change") +
    scale_fill_manual(values = c("red", "blue"), name = "Regulation", labels = c("Downregulated", "Upregulated")) +
    theme_minimal()
  
  # Save the plot
  ggsave(filename = paste0("Top_30_DEGs_", sample, ".png"), plot = plot)
}

#install.packages("plotly")
#install.packages("htmlwidgets")
#install.packages("edgeR")
#BiocManager::install("edgeR")
#BiocManager::install("ggplot2")
#BiocManager::install("plotly")
#install.packages("plotly")

# Load necessary libraries
library(DESeq2)
library(plotly)
library(edgeR)
library(htmlwidgets)
library(ggplot2)


# Function to get the top 30 upregulated and downregulated DEGs
get_top_30_up_down_DEGs <- function(res) {
  res_ordered_up <- res[order(-res$log2FoldChange, na.last = NA), ]
  res_ordered_down <- res[order(res$log2FoldChange, na.last = NA), ]
  top_30_up <- head(res_ordered_up, 30)
  top_30_down <- head(res_ordered_down, 30)
  return(list(up = top_30_up, down = top_30_down))
}

for (sample in names(results_list)) {
  res <- results_list[[sample]]
  top_30 <- get_top_30_up_down_DEGs(res)
  
  # Combine upregulated and downregulated data
  top_30_combined <- rbind(top_30$up, top_30$down)
  top_30_combined$regulation <- c(rep("Upregulated", 30), rep("Downregulated", 30))
  
  # Create an interactive bar plot for the top 30 upregulated and downregulated DEGs
  plot <- ggplot(top_30_combined, aes(x = reorder(rownames(top_30_combined), log2FoldChange), y = log2FoldChange, fill = regulation)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top 30 Upregulated and Downregulated DEGs in", sample),
         x = "Genes",
         y = "Log2 Fold Change") +
    scale_fill_manual(values = c("red", "blue"), name = "Regulation", labels = c("Downregulated", "Upregulated")) +
    theme_minimal()
  
  # Convert to plotly interactive plot
  plotly_plot <- ggplotly(plot)
  
  # Save the interactive plot as an HTML file
  htmlwidgets::saveWidget(as_widget(plotly_plot), file = paste0("Top_30_DEGs_", sample, ".html"))
}


# Function to create a volcano plot
create_volcano_plot <- function(res, sample_name) {
  # Create a data frame from the results
  res_df <- as.data.frame(res)
  
  # Add a column to categorize significant points
  res_df$Significant <- "Not Significant"
  res_df$Significant[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Significant"
  
  # Create the volcano plot
  plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.4) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    labs(title = paste("Volcano Plot for", sample_name),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal()
  
  # Save the plot
  ggsave(filename = paste0("Volcano_Plot_", sample_name, ".png"), plot = plot)
}

# Load necessary libraries
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

# Function to get top 50 DEGs
get_top_50_DEGs <- function(res) {
  res_ordered <- res[order(res$padj, na.last = NA), ]
  top_50 <- head(res_ordered, 100)
  return(rownames(top_50))
}

# Ensure metadata has row names matching the sample names
rownames(metadata) <- metadata$Samples

# Function to create and save a heatmap
create_heatmap <- function(count_data, top_genes, sample_name, metadata) {
  # Subset the count data for the top genes
  top_genes_data <- count_data[top_genes, ]
  
  # Normalize the count data
  top_genes_data <- t(apply(top_genes_data, 1, scale))
  
  # Ensure metadata rows match the columns in count data
  samples_in_data <- colnames(top_genes_data)
  annotation_data <- metadata[samples_in_data, , drop = FALSE]
  
  # Set file name
  file_name <- paste0("Heatmap_Top_50_DEGs_", sample_name, ".png")
  
  # Create the heatmap
  heatmap <- Heatmap(top_genes_data, 
                     name = "expression", 
                     row_names_gp = gpar(fontsize = 10), 
                     column_names_gp = gpar(fontsize = 10),
                     top_annotation = HeatmapAnnotation(df = annotation_data))
  
  # Save the heatmap as a PNG file
  png(filename = file_name, width = 800, height = 800)
  draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# Loop through each sample and create the heatmap
for (sample in names(results_list)) {
  res <- results_list[[sample]]
  top_50_genes <- get_top_50_DEGs(res)
  create_heatmap(count_data_subset, top_50_genes, sample, metadata)
}


# Load necessary libraries
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

# Function to filter and get top 50 DEGs
get_filtered_top_50_DEGs <- function(res) {
  # Replace NAs with 0
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  res$padj[is.na(res$padj)] <- 1
  
  res_filtered <- res[res$padj < 0.05 & (res$log2FoldChange > 1 | res$log2FoldChange < -1), ]
  res_ordered <- res_filtered[order(res_filtered$padj, na.last = NA), ]
  top_50 <- head(res_ordered, 50)
  return(rownames(top_50))
}

# Ensure metadata has row names matching the sample names
rownames(metadata) <- metadata$Samples

# Function to create and save an aggregated heatmap with sample labels
create_aggregated_heatmap_with_labels <- function(count_data, top_genes, metadata) {
  # Subset the count data for the top genes
  top_genes_data <- count_data[top_genes, ]
  
  # Normalize the count data
  top_genes_data <- t(apply(top_genes_data, 1, scale))
  
  # Ensure metadata rows match the columns in count data
  samples_in_data <- colnames(top_genes_data)
  
  # Set file name
  file_name <- "Aggregated_Heatmap_Top_50_DEGs_with_Labels.png"
  
  # Create the heatmap with sample labels
  heatmap <- Heatmap(top_genes_data, 
                     name = "expression", 
                     row_names_gp = gpar(fontsize = 10), 
                     column_names_gp = gpar(fontsize = 10),
                     column_labels = samples_in_data)  # Add sample labels
  
  # Save the heatmap as a PNG file
  png(filename = file_name, width = 2000, height = 6000)
  draw(heatmap, heatmap_legend_side = "right")
  dev.off()
}

# Get the top 50 DEGs for each sample and combine unique genes
all_top_genes <- unique(unlist(lapply(results_list, get_filtered_top_50_DEGs)))

# Create and save the aggregated heatmap with sample labels
create_aggregated_heatmap_with_labels(count_data_subset, all_top_genes, metadata)




# Load necessary libraries
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

# Function to filter and get DEGs
get_filtered_DEGs <- function(res) {
  # Replace NAs with 0
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  res$padj[is.na(res$padj)] <- 1
  
  res_filtered <- res[res$padj < 0.05 & (res$log2FoldChange > 1 | res$log2FoldChange < -1), ]
  return(res_filtered)
}

# Ensure metadata has row names matching the sample names
rownames(metadata) <- metadata$Samples

# Calculate the mean of control samples for each gene
control_samples <- rownames(metadata[metadata$Condition == "Control", ])
mean_control_values <- rowMeans(count_data_subset[, control_samples], na.rm = TRUE)

# Function to create and save individual heatmaps vs control mean
create_individual_heatmaps_vs_control_mean <- function(count_data, results_list, metadata, mean_control_values) {
  unique_conditions <- unique(metadata$Condition)
  
  for (condition in unique_conditions) {
    if (condition == "Control") next
    
    # Get samples for the current condition
    condition_samples <- rownames(metadata[metadata$Condition == condition, ])
    
    # Get filtered DEGs for the current condition
    res <- results_list[[condition]]
    filtered_res <- get_filtered_DEGs(res)
    
    # Get the top DEGs based on filtering
    top_genes <- rownames(filtered_res)
    
    # Subset the count data for the top genes and the condition samples
    sample_data <- count_data[top_genes, condition_samples, drop = FALSE]
    
    # Add the mean control values as a new column
    sample_data <- cbind(mean_control_values[top_genes], sample_data)
    colnames(sample_data)[1] <- "Mean_Control"
    
    # Remove rows with all zeros
    sample_data <- sample_data[rowSums(sample_data != 0) > 0, , drop = FALSE]
    
    # If there is no data left after removing zeros, skip this condition
    if (nrow(sample_data) == 0) {
      next
    }
    
    # Normalize the count data
    sample_data <- t(apply(sample_data, 1, scale))
    
    # Create column annotations with sample labels
    sample_labels <- c("Mean_Control", condition_samples)
    
    # Sort sample labels and data
    sorted_indices <- order(sample_labels)
    sample_labels <- sample_labels[sorted_indices]
    sample_data <- sample_data[, sorted_indices]
    
    annotation <- HeatmapAnnotation(Sample = anno_text(sample_labels))
    
    # Set file name
    file_name <- paste0("Heatmap_DEGs_", condition, "_vs_Control_Mean.png")
    
    # Determine dynamic dimensions
    num_genes <- nrow(sample_data)
    num_samples <- ncol(sample_data)
    
    # Create the heatmap
    heatmap <- Heatmap(sample_data, 
                       name = "expression", 
                       row_names_gp = gpar(fontsize = 10), 
                       column_names_gp = gpar(fontsize = 10),
                       top_annotation = annotation,
                       column_labels = sample_labels)  # Add sample labels
    
    # Save the heatmap as a PNG file with dynamic dimensions
    png(filename = file_name, width = 1000, height = 400+ (num_genes * 25))
    draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
    dev.off()
  }
}

# Create and save the individual heatmaps vs control mean
create_individual_heatmaps_vs_control_mean(count_data_subset, results_list, metadata, mean_control_values)


# Load necessary libraries
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(plotly)
library(htmlwidgets)
library(htmltools)

# Function to filter and get DEGs
get_filtered_DEGs <- function(res) {
  # Replace NAs with 0
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  res$padj[is.na(res$padj)] <- 1
  
  res_filtered <- res[res$padj < 0.05 & (res$log2FoldChange > 1 | res$log2FoldChange < -1), ]
  return(res_filtered)
}

# Ensure metadata has row names matching the sample names
rownames(metadata) <- metadata$Samples

# Calculate the mean of control samples for each gene
control_samples <- rownames(metadata[metadata$Condition == "Control", ])
mean_control_values <- rowMeans(count_data_subset[, control_samples], na.rm = TRUE)

# Function to create and save individual interactive heatmaps vs control mean
create_interactive_heatmaps_vs_control_mean <- function(count_data, results_list, metadata, mean_control_values) {
  unique_conditions <- unique(metadata$Condition)
  html_widgets <- list()
  
  for (condition in unique_conditions) {
    if (condition == "Control") next
    
    # Get samples for the current condition
    condition_samples <- rownames(metadata[metadata$Condition == condition, ])
    
    # Get filtered DEGs for the current condition
    res <- results_list[[condition]]
    filtered_res <- get_filtered_DEGs(res)
    
    # Get the top DEGs based on filtering
    top_genes <- rownames(filtered_res)
    
    # Subset the count data for the top genes and the condition samples
    sample_data <- count_data[top_genes, condition_samples, drop = FALSE]
    
    # Add the mean control values as a new column
    sample_data <- cbind(mean_control_values[top_genes], sample_data)
    colnames(sample_data)[1] <- "Mean_Control"
    
    # Remove rows with all zeros
    sample_data <- sample_data[rowSums(sample_data != 0) > 0, , drop = FALSE]
    
    # If there is no data left after removing zeros, skip this condition
    if (nrow(sample_data) == 0) {
      next
    }
    
    # Normalize the count data
    sample_data <- t(apply(sample_data, 1, scale))
    
    # Create column annotations with sample labels
    sample_labels <- c("Mean_Control", condition_samples)
    
    # Sort sample labels and data
    sorted_indices <- order(sample_labels)
    sample_labels <- sample_labels[sorted_indices]
    sample_data <- sample_data[, sorted_indices]
    
    # Create the interactive heatmap using plotly
    heatmap_plot <- plot_ly(
      z = sample_data,
      x = sample_labels,
      y = rownames(sample_data),
      type = "heatmap",
      colorscale = "Viridis"
    ) %>%
      layout(
        title = paste0("Heatmap of DEGs in ", condition, " vs Control Mean"),
        xaxis = list(title = "Samples"),
        yaxis = list(title = "Genes")
      )
    
    # Save each plotly heatmap as an HTML widget
    heatmap_widget <- htmlwidgets::saveWidget(heatmap_plot, file = paste0("Heatmap_", condition, "_vs_Control_Mean.html"), selfcontained = TRUE)
    
    html_widgets[[condition]] <- heatmap_widget
  }
  
  # Combine all individual HTML widgets into a single HTML document
  combined_html <- tagList(html_widgets)
  htmltools::save_html(combined_html, file = "Combined_Interactive_Heatmaps.html")
}

# Create and save the individual interactive heatmaps vs control mean
create_interactive_heatmaps_vs_control_mean(count_data_subset, results_list, metadata, mean_control_values)


# Install necessary packages if not already installed
# if (!requireNamespace("edgeR", quietly = TRUE)) {
#   install.packages("edgeR")
# }
# if (!requireNamespace("limma", quietly = TRUE)) {
#   install.packages("limma")
# }
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("Glimma")

# Load necessary libraries
library(edgeR)
library(limma)
library(Glimma)

# Create DGEList object from count data
dge <- DGEList(counts = count_data_subset, group = metadata$Condition)

# Calculate normalization factors
dge <- calcNormFactors(dge)

# Generate colors and pch for groups
group_colors <- rainbow(length(unique(dge$samples$group)))
names(group_colors) <- unique(dge$samples$group)
pch <- as.numeric(factor(dge$samples$group))

# Define output path
output_filepath <- "path/to/output"

# Ensure the output directory exists
if (!dir.exists(output_filepath)) {
  dir.create(output_filepath, recursive = TRUE)
}

# Assuming metadata has columns named "Samples", "Condition", and "Group"
# Modify the column names as per your actual metadata structure

# Function to generate interactive MDS plots
generate_interactive_MDS_plots <- function(dge, metadata, output_filepath) {
  # MDS plot of all genes
  glMDSPlot(cpm(dge, log = TRUE), top = Inf,
            col = group_colors[dge$samples$group],
            pch = pch[dge$samples$group], dim = c(1,2),
            groups = metadata[, c("Samples", "Condition", "Group")],
            main = 'MDS plot of all genes',
            html = 'MDS_plot_allGenes.html',
            labels = metadata$Samples,
            path = output_filepath,
            folder = 'Interactive_Plots')
  
  # MDS plot of 500 most variable genes (pairwise selection)
  glMDSPlot(cpm(dge, log = TRUE), top = 700,
            gene.selection = 'pairwise',
            col = group_colors[dge$samples$group],
            pch = pch[dge$samples$group], dim = c(1,2),
            groups = metadata[, c("Samples", "Condition", "Group")],
            main = 'MDS plot of 500 most variable genes (pairwise selection)',
            html = 'MDS_plot_500_pairwise.html',
            labels = metadata$Samples,
            path = output_filepath,
            folder = 'Interactive_Plots')
  
  # MDS plot of 500 most variable genes (common selection)
  glMDSPlot(cpm(dge, log = TRUE), top = 700,
            gene.selection = 'common',
            col = group_colors[dge$samples$group],
            pch = pch[dge$samples$group], dim = c(1,2),
            groups = metadata[, c("Samples", "Condition", "Group")],
            main = 'MDS plot of 500 most variable genes (common selection)',
            html = 'MDS_plot_500_common.html',
            labels = metadata$Samples,
            path = output_filepath,
            folder = 'Interactive_Plots')
}

# Generate the interactive MDS plots
generate_interactive_MDS_plots(dge, metadata, output_filepath)


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
# Assuming deg_counts is already defined as a data frame with columns: Sample, FoldChange1, FoldChange15, FoldChange2
##make a pie chart
bins <- c(15, 50, 150, 250, 500, 700)

# Create a factor variable with the bins
data_bins1 <- cut(deg_counts$FoldChange1, breaks = c(-Inf, bins), labels = c("0-15", "15-50", "50-150", "150-250", "250-500", "500-800"))
data_bins15 <- cut(deg_counts$FoldChange15, breaks = c(-Inf, bins), labels = c("0-15", "15-50", "50-150", "150-250", "250-500", "500-800"))
data_bins2 <- cut(deg_counts$FoldChange2, breaks = c(-Inf, bins), labels = c("0-15", "15-50", "50-150", "150-250", "250-500", "500-800"))                 
# Create a data frame with the binned data
data_frame1 <- data.frame(data = deg_counts$FoldChange1, bin = data_bins1)
data_frame15 <- data.frame(data = deg_counts$FoldChange15, bin = data_bins15)
data_frame2 <- data.frame(data = deg_counts$FoldChange2, bin = data_bins2)

# Summarize the data by bin
bin_summary1 <- data_frame1 %>%
  count(bin) %>%
  mutate(percentage = n / sum(n) * 100)
bin_summary15 <- data_frame15 %>%
  count(bin) %>%
  mutate(percentage = n / sum(n) * 100)
bin_summary2 <- data_frame2 %>%
  count(bin) %>%
  mutate(percentage = n / sum(n) * 100)
##remove the last row
bin_summary1 <- bin_summary1[-7,]
bin_summary15 <- bin_summary15[-7,]
bin_summary2 <- bin_summary2[-7,]

# Print the summary
print(bin_summary1)

# Function to create a pie chart
create_pie_chart <- function(bin_summary, title) {
  # Calculate the position for the text labels
  bin_summary <- bin_summary %>%
    arrange(desc(bin)) %>%
    mutate(ypos = cumsum(percentage) - 0.5 * percentage)
  
  # Create a pie chart
  ggplot(bin_summary, aes(x = "", y = percentage, fill = bin)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    labs(title = title, x = "", y = "", fill = "Count Bins") +
    geom_text(aes(label = paste0(round(percentage, 1), "%"), y = ypos), color = "black", size = 5) +
    theme(legend.title = element_blank()) +
    scale_fill_brewer(palette = "Set3")
}

# Create pie charts
pie_chart1 <- create_pie_chart(bin_summary1, "DEG Counts for Fold Change 1")
pie_chart15 <- create_pie_chart(bin_summary15, "DEG Counts for Fold Change 15")
pie_chart2 <- create_pie_chart(bin_summary2, "DEG Counts for Fold Change 2")

# Print the pie charts
print(pie_chart1)
print(pie_chart15)
print(pie_chart2)

ggsave("pie_chart_foldchange1.png", plot = pie_chart1, width = 8, height = 6)
ggsave("pie_chart_foldchange15.png", plot = pie_chart15, width = 8, height = 6)
ggsave("pie_chart_foldchange2.png", plot = pie_chart2, width = 8, height = 6)


library(VennDiagram)

# Function to extract gene names at different fold changes
extract_genes <- function(res, fold_change) {
  abs_log2fc <- log2(fold_change)
  genes <- rownames(res)[res$padj < 0.05 & abs(res$log2FoldChange) >= abs_log2fc]
  genes[!is.na(genes)]
}

# Initialize a data frame to store DEG details
deg_details <- data.frame(Sample = character(),
                          FoldChange1_Genes = character(),
                          FoldChange15_Genes = character(),
                          FoldChange2_Genes = character(),
                          stringsAsFactors = FALSE)

# Loop over results and extract genes
for (sample in names(results_list)) {
  res <- results_list[[sample]]
  foldchange1_genes <- paste(extract_genes(res, 1), collapse = ", ")
  foldchange15_genes <- paste(extract_genes(res, 1.5), collapse = ", ")
  foldchange2_genes <- paste(extract_genes(res, 2), collapse = ", ")
  deg_details <- rbind(deg_details, data.frame(Sample = sample,
                                               FoldChange1_Genes = foldchange1_genes,
                                               FoldChange15_Genes = foldchange15_genes,
                                               FoldChange2_Genes = foldchange2_genes))
}

print(deg_details)


# Function to create Venn diagram for each sample
create_venn_diagram <- function(sample_name, genes1, genes15, genes2) {
  venn.plot <- venn.diagram(
    x = list(
      FoldChange1 = unlist(strsplit(genes1, ", ")),
      FoldChange15 = unlist(strsplit(genes15, ", ")),
      FoldChange2 = unlist(strsplit(genes2, ", "))
    ),
    category.names = c("FoldChange1", "FoldChange15", "FoldChange2"),
    filename = paste0("VennDiagram_", sample_name, ".png"),
    output = TRUE,
    imagetype = "png",
    height = 3000,
    width = 3000,
    resolution = 100,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = c(alpha("blue", 0.5), alpha("red", 0.5), alpha("green", 0.5)),
    cex = 15,
    cat.cex = 15,
    cat.pos = 0,
    cat.dist = 0.05,
    cat.default.pos = "outer",
    cat.col = c("blue", "red", "green")
  )
}

# Create Venn diagrams for each sample
for (i in 1:nrow(deg_details)) {
  sample_name <- deg_details$Sample[i]
  genes1 <- deg_details$FoldChange1_Genes[i]
  genes15 <- deg_details$FoldChange15_Genes[i]
  genes2 <- deg_details$FoldChange2_Genes[i]
  
  create_venn_diagram(sample_name, genes1, genes15, genes2)
}









# Initialize a list to store results
results_list <- list()

# Compute results for each level (sample) vs control
levels_to_test <- setdiff(levels(metadata$Condition), "Control")
for (level in levels_to_test) {
  res <- results(dds, contrast = c("Condition", level, "Control"))
  results_list[[level]] <- res
}

# Function to identify DEGs at different fold changes
identify_DEGs <- function(res, fold_change) {
  abs_log2fc <- log2(fold_change)
  rownames(res)[res$padj < 0.05 & abs(res$log2FoldChange) >= abs_log2fc]
}
library(VennDiagram)
deg_sets <- list()

# Loop over results and identify DEGs for each fold change
for (sample in names(results_list)) {
  res <- results_list[[sample]]
  deg_sets[[sample]] <- list(
    FoldChange1 = identify_DEGs(res, 1),
    FoldChange15 = identify_DEGs(res, 15),
    FoldChange2 = identify_DEGs(res, 2)
  )
}

# Function to create Venn diagram for a sample
create_venn_diagram <- function(deg_set, sample_name) {
  venn.plot <- venn.diagram(
    x = list(
      FoldChange1 = deg_set$FoldChange1,
      FoldChange15 = deg_set$FoldChange15,
      FoldChange2 = deg_set$FoldChange2
    ),
    category.names = c("Fold Change 1", "Fold Change 15", "Fold Change 2"),
    filename = paste0("venn_", sample_name, ".png"),
    output = TRUE
  )
  return(venn.plot)
}

# Create Venn diagrams for each sample
for (sample in names(deg_sets)) {
  create_venn_diagram(deg_sets[[sample]], sample)
}
