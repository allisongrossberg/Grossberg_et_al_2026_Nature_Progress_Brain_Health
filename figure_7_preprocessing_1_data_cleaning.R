# =============================================================================
# Figure 7 – Preprocessing step 1: Proteomics data cleaning (HA lysate)
# =============================================================================
#
# Description: Loads HA lysate proteomics data, filters by detection rate,
#   replaces zeros for log2 transform, and exports long-format table for DEA.
#
# Prerequisites: figure_7_input_data/AGrossberg_COAST_study_proteomics_20240618.xlsx
#
# Output: figure_7_output_data/pro_HA_cleaned_long_data.csv
#
# Usage: Run from project root or figure_7 directory.
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}

library(tidyr)
library(tidyverse)
library(tibble)
library(readxl)

INPUT_DIR <- "figure_7_input_data"
OUTPUT_DIR <- "figure_7_output_data"
proteomics_xlsx <- "AGrossberg_COAST_study_proteomics_20240618.xlsx"
proteomics_file <- file.path(INPUT_DIR, proteomics_xlsx)
if (!file.exists(proteomics_file)) {
  stop("Proteomics data file not found. Place ", proteomics_xlsx, " in figure_7_input_data/.")
}

# Load raw data
AGrossberg_COAST_study_proteomics_20240618 <- read_excel(proteomics_file)
Pro_HA_raw_data <- AGrossberg_COAST_study_proteomics_20240618

# Select data columns
protein_sample_HA <- Pro_HA_raw_data %>% dplyr::select(2, 3, 4, 8, 135:164)

# Rename columns (remove "AGrossberg_" prefix)
colnames(protein_sample_HA) <- gsub("AGrossberg_", "", colnames(protein_sample_HA))
colnames(protein_sample_HA) <- gsub("MaxLFQ Intensity", "", colnames(protein_sample_HA))

# Remove "_HUMAN" suffix from Entry Name
protein_sample_HA$`Entry Name` <- sub("_HUMAN$", "", protein_sample_HA$`Entry Name`)

# Rename key columns
protein_sample_HA <- protein_sample_HA %>%
  dplyr::rename(Entry_Name = `Entry Name`)

# Rename Protein ID column
protein_sample_HA <- protein_sample_HA %>%
  dplyr::rename(Protein_ID = `Protein ID`)

# Convert dataframe to long format
protein_sample_HA_long <- protein_sample_HA %>%
  pivot_longer(cols = -c(Protein_ID, Entry_Name, Gene, Description),
               names_to = "Condition", 
               values_to = "MaxLFQ_Intensity")


# Step 1: Calculate the percentage of non-zero values for each protein
protein_summary <- protein_sample_HA_long %>%
  group_by(Entry_Name) %>%
  summarize(
    total_measurements = n(),
    non_zero_measurements = sum(MaxLFQ_Intensity != 0),
    non_zero_percentage = non_zero_measurements / total_measurements * 100
  )

# Step 2: Define a threshold for filtering (e.g., 50% non-zero values)
threshold <- 50

# Step 3: Get the list of proteins to keep
proteins_to_keep <- protein_summary %>%
  filter(non_zero_percentage >= threshold) %>%
  pull(Entry_Name)

# Step 4: Filter the original dataframe
protein_sample_HA_long_filtered <- protein_sample_HA_long %>%
  filter(Entry_Name %in% proteins_to_keep)

# Count zero values in filtered data
sum(protein_sample_HA_long_filtered$MaxLFQ_Intensity == 0)

# Find the smallest non-zero value
min_non_zero <- min(protein_sample_HA_long_filtered$MaxLFQ_Intensity[protein_sample_HA_long_filtered$MaxLFQ_Intensity > 0], na.rm = TRUE)
min_non_zero

# Replace zeros with half of the smallest non-zero value (for log transform)
protein_sample_HA_long_filtered_sc_log <- protein_sample_HA_long_filtered %>%
  mutate(log2_Intensity = case_when(
    MaxLFQ_Intensity > 0 ~ log2(MaxLFQ_Intensity),
    TRUE ~ log2(min_non_zero / 2)
  ))

sum(protein_sample_HA_long_filtered_sc_log$log2_Intensity == 0)

# Summary statistics
summary(protein_sample_HA_long_filtered_sc_log$log2_Intensity)
summary(protein_sample_HA_long_filtered$MaxLFQ_Intensity)

# Standard deviation
sd(protein_sample_HA_long_filtered_sc_log$log2_Intensity)

# Histogram
hist(protein_sample_HA_long_filtered_sc_log$log2_Intensity, 
     breaks = 50, 
     main = "Distribution of log2 Intensities")

# Histogram
hist(protein_sample_HA_long_filtered$MaxLFQ_Intensity, 
     breaks = 50, 
     main = "Distribution of log2 Intensities")

# Compare original vs log2 transformed
plot(protein_sample_HA_long_filtered$MaxLFQ_Intensity,
     protein_sample_HA_long_filtered_sc_log$log2_Intensity,
     xlab = "Original Intensity", 
     ylab = "Log2 Intensity",
     main = "Original vs Log2 Transformed Intensities")

qqnorm(protein_sample_HA_long_filtered_sc_log$log2_Intensity)
qqline(protein_sample_HA_long_filtered_sc_log$log2_Intensity, col = "red")

# Save cleaned proteomics data and summary dataframes to output folder
dir.create(INPUT_DIR, showWarnings = FALSE)
dir.create(OUTPUT_DIR, showWarnings = FALSE)
out_cols <- c("Entry_Name", "Protein_ID", "Condition", "log2_Intensity")
out_cols <- out_cols[out_cols %in% names(protein_sample_HA_long_filtered_sc_log)]
write.csv(protein_sample_HA_long_filtered_sc_log %>% dplyr::select(dplyr::all_of(out_cols)),
          file.path(OUTPUT_DIR, "pro_HA_cleaned_long_data.csv"), row.names = FALSE)
write.csv(as.data.frame(protein_summary), file.path(OUTPUT_DIR, "pro_HA_protein_summary.csv"), row.names = FALSE)
cat("Saved: ", file.path(OUTPUT_DIR, "pro_HA_cleaned_long_data.csv"), " (", length(out_cols), " columns)\n")
cat("Saved: ", file.path(OUTPUT_DIR, "pro_HA_protein_summary.csv"), "\n")
