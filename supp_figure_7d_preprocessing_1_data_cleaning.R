# =============================================================================
# Supplementary Figure 7d preprocessing 1 - miRNA data cleaning
# =============================================================================
#
# Description: Loads TMM-normalized mature miRNA CPM counts for the same 24
#   ADE samples used in the proteomics pipeline, reshapes to long format,
#   and joins each sample to its biological group (Cont_Exo, COV_TBI_Exo,
#   TBI_Only_Exo, COV_Only_Exo).
#
# Prerequisites: supp_figure_7_input_data/mature_normalized_CPM.xlsx
#
# Output: supp_figure_7_output_data/supp_figure_7d_cleaned_long_data.csv,
#   supp_figure_7_output_data/supp_figure_7d_group_id_mapping.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "supp_figure_7") {
  if (file.exists("supp_figure_7")) setwd("supp_figure_7")
  else if (file.exists(file.path("..", "supp_figure_7"))) setwd(file.path("..", "supp_figure_7"))
}

library(dplyr)
library(tidyr)
library(readxl)
library(tidyverse)

INPUT_DIR <- "supp_figure_7_input_data"
OUTPUT_DIR <- "supp_figure_7_output_data"
cpm_xlsx <- "mature_normalized_CPM.xlsx"
cpm_file <- file.path(INPUT_DIR, cpm_xlsx)
if (!file.exists(cpm_file)) {
  stop("miRNA CPM data file not found. Place ", cpm_xlsx, " in supp_figure_7_input_data/.")
}

mature_normalized_CPM <- read_excel(cpm_file)
mature_normalized_miRNA_R4A <- mature_normalized_CPM

mature_normalized_miRNA_R4A_long <- mature_normalized_miRNA_R4A %>%
  pivot_longer(cols = -miRNA, names_to = "ID", values_to = "miRNA_Count")

# -----------------------------------------------------------------------------
# Map each sample code (ID) to its biological group -- same 24 ADE samples
# and group assignment used in the proteomics pipeline.
# -----------------------------------------------------------------------------
group_id_mapping <- c(
  "Cont_Exo" = "8B2PC",
  "Cont_Exo" = "IPH6E",
  "Cont_Exo" = "KBPYT",
  "Cont_Exo" = "95KXN",
  "Cont_Exo" = "DTMBS",
  "Cont_Exo" = "BW3AP",
  "COV_TBI_Exo" = "A0BND",
  "COV_TBI_Exo" = "WD32N",
  "COV_TBI_Exo" = "KJHE7",
  "COV_TBI_Exo" = "XQANO",
  "COV_TBI_Exo" = "IO6KM",
  "COV_TBI_Exo" = "P89EN",
  "TBI_Only_Exo" = "PRO10",
  "TBI_Only_Exo" = "BL52P",
  "TBI_Only_Exo" = "QWMG2",
  "TBI_Only_Exo" = "84HZC",
  "TBI_Only_Exo" = "PKKZL",
  "TBI_Only_Exo" = "VXY1B",
  "COV_Only_Exo" = "B9V50",
  "COV_Only_Exo" = "ELFHW",
  "COV_Only_Exo" = "9ER2B",
  "COV_Only_Exo" = "U7DJQ",
  "COV_Only_Exo" = "74UOH",
  "COV_Only_Exo" = "UWK9R")

group_id_df <- data.frame(
  Group = names(group_id_mapping),
  ID = unname(group_id_mapping)
)

mature_normalized_miRNA_R4A_long_group <- mature_normalized_miRNA_R4A_long %>%
  left_join(group_id_df, by = "ID")

unmapped <- mature_normalized_miRNA_R4A_long_group %>% filter(is.na(Group)) %>% pull(ID) %>% unique()
if (length(unmapped) > 0) {
  stop("Sample codes with no group mapping: ", paste(unmapped, collapse = ", "))
}

dir.create(INPUT_DIR, showWarnings = FALSE)
dir.create(OUTPUT_DIR, showWarnings = FALSE)
write.csv(mature_normalized_miRNA_R4A_long_group, file.path(OUTPUT_DIR, "supp_figure_7d_cleaned_long_data.csv"), row.names = FALSE)
write.csv(group_id_df, file.path(OUTPUT_DIR, "supp_figure_7d_group_id_mapping.csv"), row.names = FALSE)
# Wide count matrix (miRNA x sample), used directly by the DEA step
write.csv(mature_normalized_miRNA_R4A, file.path(OUTPUT_DIR, "supp_figure_7d_wide_counts.csv"), row.names = FALSE)

cat("Saved: ", file.path(OUTPUT_DIR, "supp_figure_7d_cleaned_long_data.csv"), " (",
    n_distinct(mature_normalized_miRNA_R4A_long_group$miRNA), " miRNAs, ",
    n_distinct(mature_normalized_miRNA_R4A_long_group$ID), " samples)\n", sep = "")
