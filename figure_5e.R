# =============================================================================
# Figure 5e – HA TLR Inhibitor Vimentin MFI by 5-group condition (vertical bar)
# =============================================================================
#
# Description: Vertical bar plot of normalized Vimentin MFI by five broad
#   conditions (PBS Control, COVID-19/mTBI status). One-way ANOVA; Tukey
#   post-hoc.
#
# Prerequisites: figure_5_input_data/HA_TLR_Inh_*.csv
#
# Output: figure_5e.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
base_dir   <- if (dir.exists("figure_5")) "figure_5" else "."
input_dir  <- file.path(base_dir, "figure_5_input_data")
output_dir <- file.path(base_dir, "figure_5_output_data")
tif_dir    <- base_dir

library(dplyr)
library(ggplot2)
library(ggsignif)

dapi_df <- read.csv(file.path(input_dir, "HA_TLR_Inh_dapi_duplicate_df.csv"))
dapi_df <- dapi_df %>% dplyr::rename(Image_File_Clean = file)
dapi_df$Image_File_Clean <- gsub("-DAPI_Duplicate_Image_Results.csv", "", dapi_df$Image_File_Clean)
dapi_df$Image_File_Clean <- gsub("TH1020_Exo_Cont", "TH1020_Cont", dapi_df$Image_File_Clean)
mfi_df <- read.csv(file.path(input_dir, "HA_TLR_Inh_mfi_df.csv"))
mfi_df <- mfi_df %>% dplyr::rename(Image_File_Clean = Image.File)
mfi_df$Image_File_Clean <- gsub(".oir", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("duplicate_C3-", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("duplicate_2_C2-", "", mfi_df$Image_File_Clean)
mfi_df$Image_File_Clean <- gsub("TH1020_Exo_Cont", "TH1020_Cont", mfi_df$Image_File_Clean)
gfap_df <- read.csv(file.path(input_dir, "HA_TLR_Inh_gfap_duplicate_df.csv"))
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
GFAP_result_df_clean <- GFAP_result_df_1 %>%
  mutate(Condition = sub("(_[^_]+){2}$", "", Image_File_Clean), Replicate = sub(".*_([^_]+_[^_]+_[^_]+)$", "\\1", Image_File_Clean)) %>%
  mutate(Condition = gsub("_[0-9]+$", "", Condition))
Vimentin_result_df_1$Image_File_Clean <- gsub("COAST_HA_ADE_", "", Vimentin_result_df_1$Image_File_Clean)
Vimentin_result_df_clean <- Vimentin_result_df_1 %>%
  mutate(Condition = sub("(_[^_]+){2}$", "", Image_File_Clean), Replicate = sub(".*_([^_]+_[^_]+_[^_]+)$", "\\1", Image_File_Clean)) %>%
  mutate(Condition = gsub("_[0-9]+$", "", Condition))
GFAP_summary_df_1 <- GFAP_result_df_clean %>% group_by(Condition) %>% summarise(Mean = mean(Normed_GFAP_MFI), SE = sd(Normed_GFAP_MFI) / sqrt(n()))
Vimentin_summary_df_1 <- Vimentin_result_df_clean %>% group_by(Condition) %>% summarise(Mean = mean(Normed_Vimentin_MFI), SE = sd(Normed_Vimentin_MFI) / sqrt(n()))

Vimentin_summary_df_3 <- Vimentin_summary_df_1 %>%
  filter(Condition != "PEG") %>%
  filter(!grepl("_Veh$", Condition))

conditions <- c("Neg_Cont", "Exo_Cont", "Exo_TBI_Only", "Exo_COV_Only", "Exo_COV_TBI")
Vimentin_summary_df_3 <- Vimentin_summary_df_3 %>% filter(Condition %in% conditions)
Vimentin_summary_df_3$Condition <- factor(Vimentin_summary_df_3$Condition, levels = conditions)

new_labels <- c("PBS Control", "COVID-19 (-) mTBI (-) ADEs", "COVID-19 (+) mTBI (-) ADEs",
                "COVID-19 (-) mTBI (+) ADEs", "COVID-19 (+) mTBI (+) ADEs")
Vimentin_summary_df_3$Condition <- factor(Vimentin_summary_df_3$Condition,
  levels = c("Neg_Cont", "Exo_Cont", "Exo_COV_Only", "Exo_TBI_Only", "Exo_COV_TBI"), labels = new_labels)

Vimentin_result_df_filtered <- Vimentin_result_df_clean[Vimentin_result_df_clean$Condition %in% conditions, ]
# Fit ANOVA model
anova_model_Vimentin <- aov(Normed_Vimentin_MFI ~ Condition, data = Vimentin_result_df_filtered)
# Perform Tukey post-hoc tests
tukey_results_Vimentin <- TukeyHSD(anova_model_Vimentin)

condition_mapping <- c("Neg_Cont" = "PBS Control", "Exo_Cont" = "COVID-19 (-) mTBI (-) ADEs",
  "Exo_TBI_Only" = "COVID-19 (-) mTBI (+) ADEs", "Exo_COV_Only" = "COVID-19 (+) mTBI (-) ADEs",
  "Exo_COV_TBI" = "COVID-19 (+) mTBI (+) ADEs")

format_pvalue <- function(p) { if (is.na(p) || p < 0.0001) return("< 0.0001"); sprintf("%.4f", p) }

max_y <- max(Vimentin_summary_df_3$Mean + Vimentin_summary_df_3$SE, na.rm = TRUE)

p <- ggplot(Vimentin_summary_df_3, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "summary", width = 0.7, alpha = 0.95) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.7), width = 0.25, color = "black") +
  theme_minimal() +
  scale_fill_brewer(palette = "Greys", labels = new_labels) +
  theme(legend.title = element_blank(), legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid.major = element_line(color = "grey90"), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
        axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
        legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black")) +
  labs(x = "Condition", y = "Normalized Vimentin MFI") +
  scale_x_discrete(labels = c("PBS Control", "COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)",
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)"))

p <- p +
  geom_signif(y_position = c(max(Vimentin_summary_df_3$Mean) * 1.2), xmin = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_Cont"])), xmax = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_COV_TBI"])), annotation = c(format_pvalue(tukey_results_Vimentin$Condition["Exo_COV_TBI-Exo_Cont", "p adj"])), tip_length = 0, color = "black", textsize = 8, size = 1) +
  geom_signif(y_position = c(max(Vimentin_summary_df_3$Mean) * 1.4), xmin = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_Cont"])), xmax = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_TBI_Only"])), annotation = c(format_pvalue(tukey_results_Vimentin$Condition["Exo_TBI_Only-Exo_Cont", "p adj"])), tip_length = 0, color = "black", textsize = 8, size = 1) +
  geom_signif(y_position = c(max(Vimentin_summary_df_3$Mean) * 1.6), xmin = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_COV_Only"])), xmax = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_COV_TBI"])), annotation = c(format_pvalue(tukey_results_Vimentin$Condition["Exo_COV_TBI-Exo_COV_Only", "p adj"])), tip_length = 0, color = "black", textsize = 8, size = 1) +
  geom_signif(y_position = c(max(Vimentin_summary_df_3$Mean) * 1.8), xmin = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_COV_Only"])), xmax = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_TBI_Only"])), annotation = c(format_pvalue(tukey_results_Vimentin$Condition["Exo_TBI_Only-Exo_COV_Only", "p adj"])), tip_length = 0, color = "black", textsize = 8, size = 1) +
  geom_signif(y_position = c(max(Vimentin_summary_df_3$Mean) * 2), xmin = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_COV_TBI"])), xmax = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_TBI_Only"])), annotation = c(format_pvalue(tukey_results_Vimentin$Condition["Exo_TBI_Only-Exo_COV_TBI", "p adj"])), tip_length = 0, color = "black", textsize = 8, size = 1) +
  geom_signif(y_position = c(max(Vimentin_summary_df_3$Mean) * 2.2), xmin = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_COV_TBI"])), xmax = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Neg_Cont"])), annotation = c(format_pvalue(tukey_results_Vimentin$Condition["Neg_Cont-Exo_COV_TBI", "p adj"])), tip_length = 0, color = "black", textsize = 8, size = 1) +
  geom_signif(y_position = c(max(Vimentin_summary_df_3$Mean) * 2.4), xmin = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Exo_TBI_Only"])), xmax = c(which(levels(Vimentin_summary_df_3$Condition) == condition_mapping["Neg_Cont"])), annotation = c(format_pvalue(tukey_results_Vimentin$Condition["Neg_Cont-Exo_TBI_Only", "p adj"])), tip_length = 0, color = "black", textsize = 8, size = 1)

print(p)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Vimentin_summary_df_3, file.path(output_dir, "figure_5e_Vimentin_HA_TLR_Inh_5group_summary.csv"), row.names = FALSE)
write.csv(Vimentin_result_df_filtered, file.path(output_dir, "figure_5e_Vimentin_HA_TLR_Inh_5group_result.csv"), row.names = FALSE)
ggsave(file.path(tif_dir, "figure_5e.pdf"), plot = p, width = 350, height = 200, units = "mm", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_5e.pdf")), plot = p, width = 350, height = 200, units = "mm", dpi = 600, device = "tiff")
