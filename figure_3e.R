# =============================================================================
# Figure 3e – HA time series GFAP MFI bar graph
# =============================================================================
#
# Description: Bar plot of normalized GFAP MFI by time point and condition (HA
#   astrocyte-derived exosomes). Two-way ANOVA; Tukey post-hoc.
#
# Prerequisites: figure_3_input_data/HA_Time_Series_*.csv
#
# Output: figure_3e_GFAP_HA_Time_Series.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsignif)
base_dir   <- if (dir.exists("figure_3")) "figure_3" else "."
input_dir  <- file.path(base_dir, "figure_3_input_data")
output_dir <- file.path(base_dir, "figure_3_output_data")
tif_dir    <- base_dir

dapi_df <- read.csv(file.path(input_dir, "HA_Time_Series_dapi_duplicate_df.csv"))
dapi_df <- dapi_df %>% dplyr::rename(Image_File_Clean = file)
dapi_df$Image_File_Clean <- gsub("-DAPI_Duplicate_Image_Results.csv", "", dapi_df$Image_File_Clean)
mfi_df <- read.csv(file.path(input_dir, "HA_Time_Series_mfi_df.csv"))
mfi_df <- mfi_df %>% dplyr::rename(Image_File_Clean = Image.File)
mfi_df$Image_File_Clean <- gsub(".oir", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("duplicate_C3-", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("duplicate_2_C2-", "", mfi_df$Image_File_Clean)
gfap_df <- read.csv(file.path(input_dir, "HA_Time_Series_gfap_duplicate_df.csv"))
gfap_df <- gfap_df %>% dplyr::rename(Image_File_Clean = file)
gfap_df$Image_File_Clean <- gsub("-GFAP_Duplicate_Image_Results.csv", "", gfap_df$Image_File_Clean)

dapi_df_filtered <- dapi_df %>% filter(!grepl("_NfKB_", Image_File_Clean))
mfi_df_filtered <- mfi_df %>% filter(!grepl("_NfKB_", Image_File_Clean))
gfap_df_filtered <- gfap_df %>% filter(!grepl("_NfKB_", Image_File_Clean))
colnames(dapi_df_filtered)[-which(names(dapi_df_filtered) == "Image_File_Clean")] <- paste0("dapi_df_", colnames(dapi_df_filtered)[-which(names(dapi_df_filtered) == "Image_File_Clean")])
colnames(mfi_df_filtered)[-which(names(mfi_df_filtered) == "Image_File_Clean")] <- paste0("mfi_df_", colnames(mfi_df_filtered)[-which(names(mfi_df_filtered) == "Image_File_Clean")])
colnames(gfap_df_filtered)[-which(names(gfap_df_filtered) == "Image_File_Clean")] <- paste0("gfap_df_", colnames(gfap_df_filtered)[-which(names(gfap_df_filtered) == "Image_File_Clean")])
dapi_mfi_merged_df <- merge(dapi_df_filtered, mfi_df_filtered, by = "Image_File_Clean", all = TRUE)
dapi_gfap_merged_df <- merge(dapi_df_filtered, gfap_df_filtered, by = "Image_File_Clean", all = TRUE)

GFAP_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "GFAP") %>% .[!duplicated(.), ]
Vimentin_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "Vimentin") %>% .[!duplicated(.), ]
GFAP_result_df <- GFAP_df %>% mutate(Normed_GFAP_MFI = mfi_df_Mean / dapi_df_Count)
Vimentin_result_df <- Vimentin_df %>% mutate(Normed_Vimentin_MFI = mfi_df_Mean / dapi_df_Count)

GFAP_result_df_clean <- GFAP_result_df %>%
  mutate(Timepoint = str_extract(Image_File_Clean, "(?<=T)[0-9]+(_[0-9]+hr)?"),
         Timepoint = ifelse(!grepl("^T", Timepoint), paste0("T", Timepoint), Timepoint),
         Condition = ifelse(grepl("Exo_COV_TBI", Image_File_Clean), "Exo_COV_TBI", ifelse(grepl("Exo_Cont", Image_File_Clean), "Exo_Cont", "Untreated")))

GFAP_summary_df <- GFAP_result_df_clean %>% group_by(Condition, Timepoint) %>%
  summarise(Mean = mean(Normed_GFAP_MFI), SE = sd(Normed_GFAP_MFI) / sqrt(n()), .groups = "drop")

GFAP_result_df_clean_filtered <- GFAP_result_df_clean %>% filter(Timepoint != "T0")
anova_model_gfap <- aov(Normed_GFAP_MFI ~ Timepoint * Condition, data = GFAP_result_df_clean_filtered)
tukey_results_gfap <- TukeyHSD(anova_model_gfap)

format_pvalue <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p < 0.0001")
  return(paste("p =", formatC(p, format = "f", digits = 4)))
}

GFAP_summary_df_filtered <- GFAP_summary_df %>% filter(Condition != "Untreated")
new_labels <- c("COVID-19 (-) mTBI (-) ADEs", "COVID-19 (+) mTBI (+) ADEs")
max_value <- max(GFAP_summary_df_filtered$Mean, na.rm = TRUE)
GFAP_HA_TS_Plot <- ggplot(GFAP_summary_df_filtered, aes(x = Timepoint, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_brewer(palette = "Greys", labels = new_labels) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.x = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black")) +
  labs(x = "Timepoint", y = "Normalized GFAP MFI") +
  scale_x_discrete(labels = function(x) gsub("_", "-", x)) +
  coord_cartesian(ylim = c(0, max_value * 2.1)) +
  geom_signif(
    y_position = c(6, 7, 8),
    xmin = c(0.8, 0.8, 0.8),
    xmax = c(1.8, 2.8, 3.8),
    annotation = c(
      format_pvalue(tukey_results_gfap$`Timepoint:Condition`["T2_4hr:Exo_Cont-T1_2hr:Exo_Cont", "p adj"]),
      format_pvalue(tukey_results_gfap$`Timepoint:Condition`["T3_6hr:Exo_Cont-T1_2hr:Exo_Cont", "p adj"]),
      format_pvalue(tukey_results_gfap$`Timepoint:Condition`["T4_8hr:Exo_Cont-T1_2hr:Exo_Cont", "p adj"])
    ),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(9, 10, 11),
    xmin = c(1.2, 1.2, 1.2),
    xmax = c(2.2, 4.2, 5.2),
    annotation = c(
      format_pvalue(tukey_results_gfap$`Timepoint:Condition`["T2_4hr:Exo_COV_TBI-T1_2hr:Exo_COV_TBI", "p adj"]),
      format_pvalue(tukey_results_gfap$`Timepoint:Condition`["T4_8hr:Exo_COV_TBI-T1_2hr:Exo_COV_TBI", "p adj"]),
      format_pvalue(tukey_results_gfap$`Timepoint:Condition`["T5_24hr:Exo_COV_TBI-T1_2hr:Exo_COV_TBI", "p adj"])
    ),
    tip_length = 0, color = "black", textsize = 6, size = 1
  )

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(GFAP_summary_df, file.path(output_dir, "figure_3e_GFAP_HA_Time_Series_summary.csv"), row.names = FALSE)
write.csv(GFAP_result_df_clean, file.path(output_dir, "figure_3e_GFAP_HA_Time_Series_result.csv"), row.names = FALSE)
ggsave(file.path(tif_dir, "figure_3e_GFAP_HA_Time_Series.pdf"), plot = GFAP_HA_TS_Plot, width = 350, height = 200, units = "mm", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_3e_GFAP_HA_Time_Series.pdf")), plot = GFAP_HA_TS_Plot, width = 350, height = 200, units = "mm", dpi = 600, device = "tiff")
