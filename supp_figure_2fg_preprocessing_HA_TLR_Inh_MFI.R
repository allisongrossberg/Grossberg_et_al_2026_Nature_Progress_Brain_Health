# =============================================================================
# Supplementary Figure 2f/2g preprocessing – HA TLR Inhibitor experiment MFI
# =============================================================================
#
# Description: Loads and normalizes confocal imaging data from the HA 40x
#   TLR inhibitor experiment (same source data as figure_5b.R/figure_5c.R;
#   this experiment includes a PEG-only condition and a PBS negative control
#   not used in the main Figure 5 panels). Produces:
#     GFAP_result_df_clean, Vimentin_result_df_clean – per-image normalized
#       MFI with Condition labels
#     GFAP_summary_df_1, Vimentin_summary_df_1 – MFI summarised by Condition
#
# Prerequisites: ../figure_5/figure_5_input_data/HA_TLR_Inh_dapi_duplicate_df.csv,
#   HA_TLR_Inh_mfi_df.csv, HA_TLR_Inh_gfap_duplicate_df.csv
#
# Used by: supp_figure_2f.R, supp_figure_2g.R
# =============================================================================

library(dplyr)

fig5_input_dir <- "../figure_5/figure_5_input_data"

dapi_df <- read.csv(file.path(fig5_input_dir, "HA_TLR_Inh_dapi_duplicate_df.csv"))
dapi_df <- dapi_df %>% dplyr::rename(Image_File_Clean = file)
dapi_df$Image_File_Clean <- gsub("-DAPI_Duplicate_Image_Results.csv", "", dapi_df$Image_File_Clean)
dapi_df$Image_File_Clean <- gsub("TH1020_Exo_Cont", "TH1020_Cont", dapi_df$Image_File_Clean)

mfi_df <- read.csv(file.path(fig5_input_dir, "HA_TLR_Inh_mfi_df.csv"))
mfi_df <- mfi_df %>% dplyr::rename(Image_File_Clean = Image.File)
mfi_df$Image_File_Clean <- gsub(".oir", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("duplicate_C3-", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("duplicate_2_C2-", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("TH1020_Exo_Cont", "TH1020_Cont", mfi_df$Image_File_Clean)

gfap_df <- read.csv(file.path(fig5_input_dir, "HA_TLR_Inh_gfap_duplicate_df.csv"))
gfap_df <- gfap_df %>% dplyr::rename(Image_File_Clean = file)
gfap_df$Image_File_Clean <- gsub("-GFAP_Duplicate_Image_Results.csv", "", gfap_df$Image_File_Clean)
gfap_df$Image_File_Clean <- gsub("TH1020_Exo_Cont", "TH1020_Cont", gfap_df$Image_File_Clean)

colnames(dapi_df)[-which(names(dapi_df) == "Image_File_Clean")] <- paste0("dapi_df_", colnames(dapi_df)[-which(names(dapi_df) == "Image_File_Clean")])
colnames(mfi_df)[-which(names(mfi_df) == "Image_File_Clean")] <- paste0("mfi_df_", colnames(mfi_df)[-which(names(mfi_df) == "Image_File_Clean")])
colnames(gfap_df)[-which(names(gfap_df) == "Image_File_Clean")] <- paste0("gfap_df_", colnames(gfap_df)[-which(names(gfap_df) == "Image_File_Clean")])

dapi_mfi_merged_df <- merge(dapi_df, mfi_df, by = "Image_File_Clean", all = TRUE)
GFAP_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "GFAP")
Vimentin_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "Vimentin")

GFAP_result_df_1 <- GFAP_df %>% mutate(Normed_GFAP_MFI = mfi_df_Mean / dapi_df_Count)
Vimentin_result_df_1 <- Vimentin_df %>% mutate(Normed_Vimentin_MFI = mfi_df_Mean / dapi_df_Count)
GFAP_result_df_1$Image_File_Clean <- gsub("COAST_HA_ADE_", "", GFAP_result_df_1$Image_File_Clean)
Vimentin_result_df_1$Image_File_Clean <- gsub("COAST_HA_ADE_", "", Vimentin_result_df_1$Image_File_Clean)

GFAP_result_df_clean <- GFAP_result_df_1 %>%
  mutate(Condition = sub("(_[^_]+){2}$", "", Image_File_Clean), Replicate = sub(".*_([^_]+_[^_]+_[^_]+)$", "\\1", Image_File_Clean)) %>%
  mutate(Condition = gsub("_[0-9]+$", "", Condition))
Vimentin_result_df_clean <- Vimentin_result_df_1 %>%
  mutate(Condition = sub("(_[^_]+){2}$", "", Image_File_Clean), Replicate = sub(".*_([^_]+_[^_]+_[^_]+)$", "\\1", Image_File_Clean)) %>%
  mutate(Condition = gsub("_[0-9]+$", "", Condition))

GFAP_summary_df_1 <- GFAP_result_df_clean %>% group_by(Condition) %>% summarise(Mean = mean(Normed_GFAP_MFI), SE = sd(Normed_GFAP_MFI) / sqrt(n()))
Vimentin_summary_df_1 <- Vimentin_result_df_clean %>% group_by(Condition) %>% summarise(Mean = mean(Normed_Vimentin_MFI), SE = sd(Normed_Vimentin_MFI) / sqrt(n()))
