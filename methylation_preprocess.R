# ==============================================================================
# Script Name: Methylation Probe-to-Gene Matrix Processing Script
# Function: Reads 'met.gz' and 'probe_met', converts probe data to gene data 
#           (by averaging), and outputs a CSV file.
# Author: Gemini
# ==============================================================================

# 1. Check and load necessary packages
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(data.table)
library(dplyr)

# Set working directory (Ensure met.gz and probe_met are in this folder)
# setwd("C:/Your/Data/Path") 

# ==============================================================================
# Step 1: Read the probe mapping file 'probe_met'
# ==============================================================================
message(paste0("[1/5] Reading probe mapping file 'probe_met'... ", Sys.time()))

if (!file.exists("probe_met")) stop("Error: 'probe_met' file not found in the current directory!")

# Use fread for fast reading
probe_map <- fread("probe_met", header = FALSE, sep = "\t")

# Rename columns (Assuming col1 is ProbeID, col2 is GeneSymbol)
names(probe_map) <- c("ProbeID", "GeneSymbol")

# Data Cleaning: Remove rows with empty gene symbols
probe_map <- probe_map[GeneSymbol != "" & !is.na(GeneSymbol)]
# Remove duplicate mappings to ensure uniqueness
probe_map <- unique(probe_map)

message(paste("      Mapping file loaded. Total probes:", nrow(probe_map)))

# ==============================================================================
# Step 2: Read the large methylation data file 'met.gz'
# ==============================================================================
message(paste0("[2/5] Reading methylation data file 'met.gz' (Large file, please wait)... ", Sys.time()))

if (!file.exists("met.gz")) stop("Error: 'met.gz' file not found in the current directory!")

# data.table::fread is the fastest way to read large files
met_data <- fread("met.gz", header = TRUE, sep = "\t")

# Rename the first column to "ProbeID" for merging (assuming it's the probe column)
first_col_name <- colnames(met_data)[1]
setnames(met_data, first_col_name, "ProbeID")

message(paste("      Methylation data loaded. Rows:", nrow(met_data), " Columns:", ncol(met_data)))

# ==============================================================================
# Step 3: Merge Data (Probe -> Gene)
# ==============================================================================
message(paste0("[3/5] Matching probes to genes... ", Sys.time()))

# Use merge (inner join) to keep only probes that map to a gene
# This automatically filters out probes without gene annotation
merged_data <- merge(probe_map, met_data, by = "ProbeID", all = FALSE)

message(paste("      Matching complete. Valid probe rows:", nrow(merged_data)))

# ==============================================================================
# Step 4: Aggregate Data (Calculate Mean by Gene)
# ==============================================================================
message(paste0("[4/5] Aggregating data by gene (calculating mean)... This may take a while. ", Sys.time()))

# Remove the ProbeID column as we are aggregating by GeneSymbol
merged_data[, ProbeID := NULL]

# Use data.table's fast aggregation
# .SD represents all data columns; by = GeneSymbol groups by gene
# mean(..., na.rm=TRUE) calculates the average ignoring missing values
gene_matrix <- merged_data[, lapply(.SD, mean, na.rm = TRUE), by = GeneSymbol]

message(paste("      Aggregation complete! Total unique genes:", nrow(gene_matrix)))

# ==============================================================================
# Step 5: Output Result
# ==============================================================================
output_file <- "Methylation_Gene_Matrix.csv"
message(paste0("[5/5] Writing output file: ", output_file, " ... ", Sys.time()))

fwrite(gene_matrix, output_file)

message("========================================================")
message("[Success] All tasks completed!")
message(paste("File saved to:", file.path(getwd(), output_file)))
message("========================================================")