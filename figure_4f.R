# =============================================================================
# Figure 4f – GFAP / Vimentin MFI correlogram
# =============================================================================
#
# Description: Correlation matrix of GFAP MFI, Vimentin MFI, and clinical
#   variables (Pearson; FDR-adjusted p-values). Exported as TIFF.
#
# Prerequisites: figure_4_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv,
#   figure_4_input_data/RA_Ind_ADE_*.csv
#
# Output: figure_4f_gfap_vim_mfi_corr_plot.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
library(dplyr)
library(tidyverse)
base_dir   <- if (dir.exists("figure_4")) "figure_4" else "."
input_dir  <- file.path(base_dir, "figure_4_input_data")
output_dir <- file.path(base_dir, "figure_4_output_data")
tif_dir    <- base_dir

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(file.path(input_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"), stringsAsFactors = FALSE)
dapi_df_RA_Ind_ADEs <- read.csv(file.path(input_dir, "RA_Ind_ADE_dapi_duplicate_df.csv"))
dapi_df_RA_Ind_ADEs <- dapi_df_RA_Ind_ADEs %>% dplyr::rename(Image_File_Clean = file)
dapi_df_RA_Ind_ADEs$Image_File_Clean <- gsub("-DAPI_Duplicate_Image_Results.csv", "", dapi_df_RA_Ind_ADEs$Image_File_Clean)
dapi_df_RA_Ind_ADEs$Image_File_Clean <- gsub("XQAN0", "XQANO", dapi_df_RA_Ind_ADEs$Image_File_Clean)
mfi_df_RA_Ind_ADEs <- read.csv(file.path(input_dir, "RA_Ind_ADE_mfi_df.csv"))
mfi_df_RA_Ind_ADEs <- mfi_df_RA_Ind_ADEs %>% dplyr::rename(Image_File_Clean = Image.File)
mfi_df_RA_Ind_ADEs$Image_File_Clean <- gsub(".oir", "", mfi_df_RA_Ind_ADEs$Image_File_Clean)
mfi_df_RA_Ind_ADEs$Image_File_Clean <- gsub("duplicate_C3-", "", mfi_df_RA_Ind_ADEs$Image_File_Clean)
mfi_df_RA_Ind_ADEs$Image_File_Clean <- gsub("duplicate_2_C2-", "", mfi_df_RA_Ind_ADEs$Image_File_Clean)
mfi_df_RA_Ind_ADEs$Image_File_Clean <- gsub("XQAN0", "XQANO", mfi_df_RA_Ind_ADEs$Image_File_Clean)
gfap_df_RA_Ind_ADEs <- read.csv(file.path(input_dir, "RA_Ind_ADE_gfap_duplicate_df.csv"))
gfap_df_RA_Ind_ADEs <- gfap_df_RA_Ind_ADEs %>% dplyr::rename(Image_File_Clean = file)
gfap_df_RA_Ind_ADEs$Image_File_Clean <- gsub("-GFAP_Duplicate_Image_Results.csv", "", gfap_df_RA_Ind_ADEs$Image_File_Clean)
gfap_df_RA_Ind_ADEs$Image_File_Clean <- gsub("XQAN0", "XQANO", gfap_df_RA_Ind_ADEs$Image_File_Clean)
colnames(dapi_df_RA_Ind_ADEs)[-which(names(dapi_df_RA_Ind_ADEs) == "Image_File_Clean")] <- paste0("dapi_df_", colnames(dapi_df_RA_Ind_ADEs)[-which(names(dapi_df_RA_Ind_ADEs) == "Image_File_Clean")])
colnames(mfi_df_RA_Ind_ADEs)[-which(names(mfi_df_RA_Ind_ADEs) == "Image_File_Clean")] <- paste0("mfi_df_", colnames(mfi_df_RA_Ind_ADEs)[-which(names(mfi_df_RA_Ind_ADEs) == "Image_File_Clean")])
dapi_mfi_merged_df <- merge(dapi_df_RA_Ind_ADEs, mfi_df_RA_Ind_ADEs, by = "Image_File_Clean", all = TRUE)
GFAP_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "GFAP")
Vimentin_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "Vimentin")
GFAP_result_df <- GFAP_df %>% mutate(Normed_GFAP_MFI = mfi_df_Mean / dapi_df_Count)
Vimentin_result_df <- Vimentin_df %>% mutate(Normed_Vimentin_MFI = mfi_df_Mean / dapi_df_Count)
GFAP_result_df <- GFAP_result_df %>% filter(!grepl("TNFa|TGFb", Image_File_Clean))
Vimentin_result_df <- Vimentin_result_df %>% filter(!grepl("TNFa|TGFb", Image_File_Clean))
id_lists <- list(
  list(ids = c("8B2PC", "IPH6E", "KBPYT", "95KXN", "DTMBS", "BW3AP"), group = "COVID-19 (-) mTBI (-) ADEs"),
  list(ids = c("A0BND", "WD32N", "KJHE7", "XQANO", "I06KM", "P89EN"), group = "COVID-19 (+) mTBI (+) ADEs"),
  list(ids = c("PRO10", "BL52P", "QWMG2", "84HZC", "PKXZL", "VXY1B"), group = "COVID-19 (-) mTBI (+) ADEs"),
  list(ids = c("B9V50", "ELFHW", "9ER2B", "U7DJQ", "74UOH", "UWK9R"), group = "COVID-19 (+) mTBI (-) ADEs"),
  list(ids = c("LPS"), group = "LPS"),
  list(ids = c("Neg_Cont"), group = "PBS Control"))
GFAP_result_df$Group <- ""
Vimentin_result_df$Group <- ""
for (id_list in id_lists) { GFAP_result_df$Group[grepl(paste(id_list$ids, collapse = "|"), GFAP_result_df$Image_File_Clean)] <- id_list$group }
for (id_list in id_lists) { Vimentin_result_df$Group[grepl(paste(id_list$ids, collapse = "|"), Vimentin_result_df$Image_File_Clean)] <- id_list$group }
GFAP_result_df <- GFAP_result_df %>% mutate(ID = str_split(Image_File_Clean, "_") %>% sapply(function(x) ifelse(length(x) >= 4, x[4], NA)))
Vimentin_result_df <- Vimentin_result_df %>% mutate(ID = str_split(Image_File_Clean, "_") %>% sapply(function(x) ifelse(length(x) >= 4, x[4], NA)))
GFAP_summary_df <- GFAP_result_df %>% group_by(ID, Group) %>% summarise(Mean_GFAP = mean(Normed_GFAP_MFI), SE_GFAP = sd(Normed_GFAP_MFI) / sqrt(n()), .groups = "drop")
Vimentin_summary_df <- Vimentin_result_df %>% group_by(ID, Group) %>% summarise(Mean_vim = mean(Normed_Vimentin_MFI), SE_vim = sd(Normed_Vimentin_MFI) / sqrt(n()), .groups = "drop")

# -----------------------------------------------------------------------------
# Correlations
# -----------------------------------------------------------------------------
merged_gfap_vim_mfi_df <- left_join(GFAP_summary_df, Vimentin_summary_df, by = "ID") %>%
  dplyr::select(ID, Mean_GFAP, SE_GFAP, Mean_vim, SE_vim) %>%
  rename(participant_id = ID)

merged_gfap_vim_mfi_df <- merged_gfap_vim_mfi_df %>%
  filter(!(participant_id %in% c("LPS", "Neg")))

merged_gfap_vim_mfi_df <- merged_gfap_vim_mfi_df %>%
  mutate(participant_id = case_when(
    participant_id == "74UOH" ~ "74U0H",
    participant_id == "I06KM" ~ "IO6KM",
    TRUE ~ participant_id
  ))

COAST_Study_Data_Clean_Age_Groups_add_dates$average_years_since_tbi

merged_gfap_vim_mfi_df_sym <- left_join(merged_gfap_vim_mfi_df, COAST_Study_Data_Clean_Age_Groups_add_dates, by = "participant_id") %>%
  dplyr::select(Mean_GFAP, SE_GFAP, Mean_vim, SE_vim, qq_group, age_years, 
         recode_overall_neuro_psych_sev_score, recode_overall_neuro_psych_freq_score, 
         qq_tbi_num, covid_pos_test_num, average_years_since_tbi, average_years_since_covid)

merged_gfap_vim_mfi_df_sym <- merged_gfap_vim_mfi_df_sym %>%
  mutate(across(c(average_years_since_tbi, average_years_since_covid), ~ifelse(is.nan(.), NA, .)))

merged_gfap_vim_mfi_df_sym$Total_mTBI_COVID_Incidences <- merged_gfap_vim_mfi_df_sym$average_years_since_tbi + merged_gfap_vim_mfi_df_sym$average_years_since_covid

# Load required libraries
library(dplyr)
library(tidyr)
library(purrr)
library(corrplot)
library(psych)
library(DT)

# Define numeric columns for analysis
numeric_cols <- c("Mean_GFAP", "Mean_vim", "age_years", 
                  "recode_overall_neuro_psych_sev_score", "recode_overall_neuro_psych_freq_score", 
                  "qq_tbi_num", "covid_pos_test_num", "average_years_since_tbi", "average_years_since_covid", 
                  "Total_mTBI_COVID_Incidences")

# 1. Perform correlations on the entire dataset
cor_matrix <- cor(merged_gfap_vim_mfi_df_sym[numeric_cols], method = "pearson")
cor_test <- psych::corr.test(merged_gfap_vim_mfi_df_sym[numeric_cols], method = "pearson", adjust = "fdr")

# 2. Create a table with correlation coefficients and p-values
result_table <- data.frame(
  var1 = rownames(cor_matrix)[row(cor_matrix)[upper.tri(cor_matrix)]],
  var2 = colnames(cor_matrix)[col(cor_matrix)[upper.tri(cor_matrix)]],
  correlation = cor_matrix[upper.tri(cor_matrix)],
  p_value = cor_test$p[upper.tri(cor_test$p)],
  p_adjusted = cor_test$p[upper.tri(cor_test$p)]  # p-values are already adjusted
) %>%
  mutate(
    correlation = round(correlation, 3),
    p_value = round(p_value, 3),
    p_adjusted = round(p_adjusted, 3)
  ) %>%
  arrange(var1, var2)

# 3. Display the result table in an interactive viewer
datatable(result_table, 
          options = list(pageLength = 25, 
                         scrollX = TRUE, 
                         scrollY = "400px", 
                         dom = 'Bfrtip',
                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
          filter = 'top',
          caption = "Overall Correlation Results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(merged_gfap_vim_mfi_df_sym, file.path(output_dir, "figure_4f_merged_GFAP_Vim_MFI_symptom_data.csv"), row.names = FALSE)
write.csv(result_table, file.path(output_dir, "figure_4f_correlation_table.csv"), row.names = FALSE)

# Define numeric columns and their corresponding labels
numeric_cols <- c("Mean_GFAP", "Mean_vim",   "age_years", 
                  "recode_overall_neuro_psych_sev_score", "recode_overall_neuro_psych_freq_score", 
                  "qq_tbi_num", "covid_pos_test_num", 
                  "average_years_since_tbi", "average_years_since_covid")
custom_labels <- c("Mean GFAP MFI", "Mean Vimentin MFI", "Age",
                   "Neuro/Psych Severity Score", "Neuro/Psych Frequency Score", 
                   "Number of mTBIs", "Number of COVID-19 Infections", 
                   "Average Years Since Last mTBI", "Average Years Since Last COVID-19 Infection")

# Create a new dataframe with renamed columns
df_renamed <- merged_gfap_vim_mfi_df_sym[numeric_cols]
colnames(df_renamed) <- custom_labels

# Perform correlations on the renamed dataset
cor_matrix <- cor(df_renamed, method = "pearson", use = "pairwise.complete.obs")
cor_test <- psych::corr.test(df_renamed, method = "pearson", adjust = "fdr", use = "pairwise.complete.obs")

draw_corr_plot <- function() {
  par(font = 2)
  corrplot(cor_matrix,
           method = "circle",
           type = "upper",
           order = "hclust",
           addCoef.col = "black",
           tl.col = "black",
           tl.srt = 45,
           diag = FALSE,
           tl.cex = 1,
           font = 2,
           number.cex = 1,
           col = colorRampPalette(c("#2C5F9E", "white", "#C1440E"))(200),
           addgrid.col = "grey50",
           cl.pos = "r",
           cl.ratio = 0.15,
           cl.align = "r",
           cl.cex = 1.3,
           cl.offset = -0.8)
}

pdf_path <- file.path(tif_dir, "figure_4f_gfap_vim_mfi_corr_plot.pdf")
cairo_pdf(pdf_path, width = 300 / 25.4, height = 240 / 25.4)
draw_corr_plot()
dev.off()

tif_path <- sub("\\.pdf$", ".tif", pdf_path)
tiff(tif_path, width = 300 / 25.4, height = 240 / 25.4, units = "in", res = 600, compression = "lzw")
draw_corr_plot()
dev.off()
