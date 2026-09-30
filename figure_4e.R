# =============================================================================
# Figure 4e – RA individual ADE Vimentin MFI bar graph (group averages)
# =============================================================================
#
# Description: Bar plot of mean (SE) Vimentin MFI by group for RA individual
#   ADEs. One-way ANOVA; Tukey post-hoc.
#
# Prerequisites: figure_4_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv,
#   figure_4_input_data/RA_Ind_ADE_*.csv
#
# Output: figure_4e_RA_Ind_Vimentin_grouped.pdf
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

Vimentin_summary_df_2 <- Vimentin_result_df %>%
  group_by(Group) %>%
  summarise(Mean = mean(Normed_Vimentin_MFI),
            SE = sd(Normed_Vimentin_MFI) / sqrt(n()))


# Fit ANOVA model
anova_model <- aov(Normed_Vimentin_MFI ~ Group, data = Vimentin_result_df)
summary(anova_model)
# Perform Tukey post-hoc tests
tukey_results <- TukeyHSD(anova_model)
# Print Tukey post-hoc results
print(tukey_results)

Vimentin_summary_df_2$Group <- factor(Vimentin_summary_df_2$Group, levels = c("PBS Control", "LPS", "COVID-19 (-) mTBI (-) ADEs", "COVID-19 (+) mTBI (-) ADEs", 
                                                                              "COVID-19 (-) mTBI (+) ADEs", "COVID-19 (+) mTBI (+) ADEs"))

max_y <- max(Vimentin_summary_df_2$Mean + Vimentin_summary_df_2$SE, na.rm = TRUE)

# Create a stacked bar graph
RA_Ind_group_Vimentin_avg_plot <- ggplot(Vimentin_summary_df_2, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha=0.95) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.7), width = 0.25, color="black") +
  theme_minimal()+
  scale_fill_brewer(palette = "Greys") +  # Add labels here
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
       y = "Normalized Vimentin MFI",
       caption = "F(5, 1553) = 39.08, p < 0.0010")+
  scale_x_discrete(labels = c("PBS Control", "LPS", "COVID-19 (-)\nmTBI (-) ADEs", "COVID-19 (+)\nmTBI (-) ADEs", 
                              "COVID-19 (-)\nmTBI (+) ADEs", "COVID-19 (+)\nmTBI (+) ADEs")) +
  geom_signif(
    y_position = c(max_y * 1.2),
    xmin = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (-) mTBI (+) ADEs")),
    xmax = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (-) mTBI (-) ADEs")),
    annotation = c(format_pvalue(tukey_results$Group["COVID-19 (-) mTBI (+) ADEs-COVID-19 (-) mTBI (-) ADEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max_y * 1.5),
    xmin = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (+) mTBI (+) ADEs")),
    xmax = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (-) mTBI (-) ADEs")),
    annotation = c(format_pvalue(tukey_results$Group["COVID-19 (+) mTBI (+) ADEs-COVID-19 (-) mTBI (-) ADEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max_y * 1.8),
    xmin = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (+) mTBI (-) ADEs")),
    xmax = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (-) mTBI (+) ADEs")),
    annotation = c(format_pvalue(tukey_results$Group["COVID-19 (+) mTBI (-) ADEs-COVID-19 (-) mTBI (+) ADEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max_y * 2.1),
    xmin = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (+) mTBI (+) ADEs")),
    xmax = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (-) mTBI (+) ADEs")),
    annotation = c(format_pvalue(tukey_results$Group["COVID-19 (+) mTBI (+) ADEs-COVID-19 (-) mTBI (+) ADEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max_y * 2.4),
    xmin = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (+) mTBI (+) ADEs")),
    xmax = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (+) mTBI (-) ADEs")),
    annotation = c(format_pvalue(tukey_results$Group["COVID-19 (+) mTBI (+) ADEs-COVID-19 (+) mTBI (-) ADEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max_y * 2.7),
    xmin = c(which(levels(Vimentin_summary_df_2$Group) == "LPS")),
    xmax = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (+) mTBI (+) ADEs")),
    annotation = c(format_pvalue(tukey_results$Group["LPS-COVID-19 (+) mTBI (+) ADEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  geom_signif(
    y_position = c(max_y * 3),
    xmin = c(which(levels(Vimentin_summary_df_2$Group) == "PBS Control")),
    xmax = c(which(levels(Vimentin_summary_df_2$Group) == "COVID-19 (+) mTBI (+) ADEs")),
    annotation = c(format_pvalue(tukey_results$Group["PBS Control-COVID-19 (+) mTBI (+) ADEs", "p adj"])),
    tip_length = 0, color = "black", textsize = 6, size = 1
  ) +
  coord_cartesian(clip = "off", ylim = c(0, max_y * 3))


dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Vimentin_summary_df_2, file.path(output_dir, "figure_4e_RA_Ind_Vimentin_group_summary.csv"), row.names = FALSE)
write.csv(Vimentin_result_df, file.path(output_dir, "figure_4e_RA_Ind_Vimentin_result.csv"), row.names = FALSE)
ggsave(file.path(tif_dir, "figure_4e_RA_Ind_Vimentin_grouped.pdf"), plot = RA_Ind_group_Vimentin_avg_plot, width = 350, height = 200, units = "mm", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_4e_RA_Ind_Vimentin_grouped.pdf")), plot = RA_Ind_group_Vimentin_avg_plot, width = 350, height = 200, units = "mm", dpi = 600, device = "tiff")
