# =============================================================================
# Figure 4h – RA pooled MDE GFAP MFI bar graph (group averages)
# =============================================================================
#
# Description: Bar plot of mean (SE) GFAP MFI by condition (RA pooled MDEs).
#   One-way ANOVA; Tukey post-hoc.
#
# Prerequisites: figure_4_input_data/RA_Pooled_MDEs_*.csv
#
# Output: figure_4h_GFAP_RA_Pooled_MDEs.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
library(dplyr)
library(tidyverse)
library(ggsignif)
base_dir   <- if (dir.exists("figure_4")) "figure_4" else "."
input_dir  <- file.path(base_dir, "figure_4_input_data")
output_dir <- file.path(base_dir, "figure_4_output_data")
tif_dir    <- base_dir

format_pvalue <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p < 0.0001")
  return(paste("p =", formatC(p, format = "f", digits = 4)))
}
dapi_df_RA_pooled_MDEs_df <- read.csv(file.path(input_dir, "RA_Pooled_MDEs_dapi_duplicate_df.csv"))
dapi_df_RA_pooled_MDEs_df <- dapi_df_RA_pooled_MDEs_df %>% dplyr::rename(Image_File_Clean = file)
dapi_df_RA_pooled_MDEs_df$Image_File_Clean <- gsub("-DAPI_Duplicate_Image_Results.csv", "", dapi_df_RA_pooled_MDEs_df$Image_File_Clean)
mfi_df_RA_pooled_MDEs_df <- read.csv(file.path(input_dir, "RA_Pooled_MDEs_mfi_df.csv"))
mfi_df_RA_pooled_MDEs_df <- mfi_df_RA_pooled_MDEs_df %>% dplyr::rename(Image_File_Clean = Image.File)
mfi_df_RA_pooled_MDEs_df$Image_File_Clean <- gsub(".oir", "", mfi_df_RA_pooled_MDEs_df$Image_File_Clean)
mfi_df_RA_pooled_MDEs_df$Image_File_Clean <- gsub("duplicate_C3-", "", mfi_df_RA_pooled_MDEs_df$Image_File_Clean)
mfi_df_RA_pooled_MDEs_df$Image_File_Clean <- gsub("duplicate_2_C2-", "", mfi_df_RA_pooled_MDEs_df$Image_File_Clean)
gfap_df_RA_pooled_MDEs_df <- read.csv(file.path(input_dir, "RA_Pooled_MDEs_gfap_duplicate_df.csv"))
gfap_df_RA_pooled_MDEs_df <- gfap_df_RA_pooled_MDEs_df %>% dplyr::rename(Image_File_Clean = file)
gfap_df_RA_pooled_MDEs_df$Image_File_Clean <- gsub("-GFAP_Duplicate_Image_Results.csv", "", gfap_df_RA_pooled_MDEs_df$Image_File_Clean)
colnames(dapi_df_RA_pooled_MDEs_df)[-which(names(dapi_df_RA_pooled_MDEs_df) == "Image_File_Clean")] <- paste0("dapi_df_", colnames(dapi_df_RA_pooled_MDEs_df)[-which(names(dapi_df_RA_pooled_MDEs_df) == "Image_File_Clean")])
colnames(mfi_df_RA_pooled_MDEs_df)[-which(names(mfi_df_RA_pooled_MDEs_df) == "Image_File_Clean")] <- paste0("mfi_df_", colnames(mfi_df_RA_pooled_MDEs_df)[-which(names(mfi_df_RA_pooled_MDEs_df) == "Image_File_Clean")])
colnames(gfap_df_RA_pooled_MDEs_df)[-which(names(gfap_df_RA_pooled_MDEs_df) == "Image_File_Clean")] <- paste0("gfap_df_", colnames(gfap_df_RA_pooled_MDEs_df)[-which(names(gfap_df_RA_pooled_MDEs_df) == "Image_File_Clean")])
dapi_mfi_merged_df <- merge(dapi_df_RA_pooled_MDEs_df, mfi_df_RA_pooled_MDEs_df, by = "Image_File_Clean", all = TRUE)
dapi_gfap_merged_df <- merge(dapi_df_RA_pooled_MDEs_df, gfap_df_RA_pooled_MDEs_df, by = "Image_File_Clean", all = TRUE)
if ("dapi_df_Area" %in% names(dapi_gfap_merged_df)) dapi_gfap_merged_df$dapi_df_Area <- as.numeric(dapi_gfap_merged_df$dapi_df_Area)
dapi_gfap_merged_df$dapi_df_Count <- as.numeric(dapi_gfap_merged_df$dapi_df_Count)
GFAP_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "GFAP")
Vimentin_df <- dapi_mfi_merged_df %>% filter(mfi_df_Marker == "Vimentin")
GFAP_result_df <- GFAP_df %>% mutate(Normed_GFAP_MFI = mfi_df_Mean / dapi_df_Count)
Vimentin_result_df <- Vimentin_df %>% mutate(Normed_Vimentin_MFI = mfi_df_Mean / dapi_df_Count)
GFAP_result_df <- GFAP_result_df %>% filter(!str_detect(Image_File_Clean, "TNFa_IL1b") & !str_detect(Image_File_Clean, "TGFb"))
Vimentin_result_df <- Vimentin_result_df %>% filter(!str_detect(Image_File_Clean, "TNFa_IL1b") & !str_detect(Image_File_Clean, "TGFb"))
GFAP_result_df_clean <- GFAP_result_df %>%
  mutate(Condition = ifelse(grepl("Exo_COV_TBI", Image_File_Clean), "COVID-19 (+) mTBI (+) MDEs",
    ifelse(grepl("Exo_COV_Only", Image_File_Clean), "COVID-19 (+) mTBI (-) MDEs",
    ifelse(grepl("Exo_TBI_Only", Image_File_Clean), "COVID-19 (-) mTBI (+) MDEs",
    ifelse(grepl("Exo_Cont", Image_File_Clean), "COVID-19 (-) mTBI (-) MDEs",
    ifelse(grepl("Neg_Cont", Image_File_Clean), "PBS Control", "LPS"))))))
Vimentin_result_df_clean <- Vimentin_result_df %>%
  mutate(Condition = ifelse(grepl("Exo_COV_TBI", Image_File_Clean), "COVID-19 (+) mTBI (+) MDEs",
    ifelse(grepl("Exo_COV_Only", Image_File_Clean), "COVID-19 (+) mTBI (-) MDEs",
    ifelse(grepl("Exo_TBI_Only", Image_File_Clean), "COVID-19 (-) mTBI (+) MDEs",
    ifelse(grepl("Exo_Cont", Image_File_Clean), "COVID-19 (-) mTBI (-) MDEs",
    ifelse(grepl("Neg_Cont", Image_File_Clean), "PBS Control", "LPS"))))))
GFAP_summary_df <- GFAP_result_df_clean %>% group_by(Condition) %>% summarise(Mean = mean(Normed_GFAP_MFI), SE = sd(Normed_GFAP_MFI) / sqrt(n()), .groups = "drop")
Vimentin_summary_df <- Vimentin_result_df_clean %>% group_by(Condition) %>% summarise(Mean = mean(Normed_Vimentin_MFI), SE = sd(Normed_Vimentin_MFI) / sqrt(n()), .groups = "drop")

library(ggthemes)
library(ggsignif)

# Fit ANOVA model
anova_model_gfap <- aov(Normed_GFAP_MFI ~ Condition, data = GFAP_result_df_clean)
summary(anova_model_gfap)
# Perform Tukey post-hoc tests
tukey_results_gfap <- TukeyHSD(anova_model_gfap)
# Print Tukey post-hoc results
print(tukey_results_gfap)

max_y <- max(GFAP_summary_df$Mean + GFAP_summary_df$SE, na.rm = TRUE)

GFAP_summary_df$Condition <- factor(GFAP_summary_df$Condition, levels = c(
  "PBS Control", "COVID-19 (-) mTBI (-) MDEs", "COVID-19 (+) mTBI (-) MDEs", 
  "COVID-19 (-) mTBI (+) MDEs", "COVID-19 (+) mTBI (+) MDEs", "LPS"
))
# Create a stacked bar graph
GFAP_MDE_RA_plot <- ggplot(GFAP_summary_df, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha=0.95) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha=0.95) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.7), width = 0.25, color="black") +
  theme_minimal()+
  scale_fill_brewer(palette = "Greys") +  # Add labels here
  scale_x_discrete(labels = c("PBS Control", "COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)", 
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)", "LPS")) +
  theme(legend.title=element_blank(),legend.position="none", 
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.text.x = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black"), 
        axis.text.y = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black"), 
        legend.text = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black"),
        plot.caption = element_text(size = 18, margin = margin(t = 20), family = "sans", face = "bold"))+
  labs(x = "Condition",
       y = "Normalized GFAP MFI",
       caption = "F(5, 354) = 16.42, p < 0.0010") +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 1.25),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (-) mTBI (-) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (+) mTBI (+) MDEs")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["COVID-19 (+) mTBI (+) MDEs-COVID-19 (-) mTBI (-) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 1.4),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (-) mTBI (-) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "LPS")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["LPS-COVID-19 (-) mTBI (-) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 1.55),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (-) mTBI (+) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (+) mTBI (+) MDEs")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["COVID-19 (+) mTBI (+) MDEs-COVID-19 (-) mTBI (+) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 1.7),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (-) mTBI (+) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "LPS")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["LPS-COVID-19 (-) mTBI (+) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 1.85),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (+) mTBI (-) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (+) mTBI (+) MDEs")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["COVID-19 (+) mTBI (+) MDEs-COVID-19 (+) mTBI (-) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 2.0),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (+) mTBI (-) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "LPS")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["LPS-COVID-19 (+) mTBI (-) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 2.15),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (+) mTBI (+) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "LPS")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["LPS-COVID-19 (+) mTBI (+) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 2.3),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "COVID-19 (+) mTBI (+) MDEs")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "PBS Control")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["PBS Control-COVID-19 (+) mTBI (+) MDEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max(GFAP_summary_df$Mean) * 2.45),
    xmin = c(which(levels(GFAP_summary_df$Condition) == "LPS")),
    xmax = c(which(levels(GFAP_summary_df$Condition) == "PBS Control")),
    annotation = c(format_pvalue(tukey_results_gfap$Condition["PBS Control-LPS", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 0.5
  ) +
  coord_cartesian(clip = "off", ylim = c(0, max(GFAP_summary_df$Mean) * 2.5))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(GFAP_summary_df, file.path(output_dir, "figure_4h_GFAP_RA_Pooled_MDEs_summary.csv"), row.names = FALSE)
write.csv(GFAP_result_df_clean, file.path(output_dir, "figure_4h_GFAP_RA_Pooled_MDEs_result.csv"), row.names = FALSE)
ggsave(file.path(tif_dir, "figure_4h_GFAP_RA_Pooled_MDEs.pdf"), plot = GFAP_MDE_RA_plot, width = 350, height = 200, units = "mm", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_4h_GFAP_RA_Pooled_MDEs.pdf")), plot = GFAP_MDE_RA_plot, width = 350, height = 200, units = "mm", dpi = 600, device = "tiff")
