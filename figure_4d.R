# =============================================================================
# Figure 4d – RA individual ADE Vimentin MFI bar graph (per-participant)
# =============================================================================
#
# Description: Per-participant bar plot of Vimentin MFI for RA individual ADEs,
#   faceted by group.
#
# Prerequisites: figure_4_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv,
#   figure_4_input_data/RA_Ind_ADE_*.csv
#
# Output: figure_4d_RA_Ind_Vimentin.pdf
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
# Vimentin per-participant bar plot
# -----------------------------------------------------------------------------
library(ggthemes)
library(forcats)

# First, we need to reorder the Group factor levels based on the specified order

mapping <- c(
  "8B2PC" = 1, "IPH6E" = 6, "KBPYT" = 5, "95KXN" = 4, "DTMBS" = 2, "BW3AP" = 3,
  "A0BND" = 1, "WD32N" = 2, "KJHE7" = 3, "XQANO" = 4, "I06KM" = 5, "P89EN" = 6,
  "PRO10" = 1, "BL52P" = 2, "QWMG2" = 3, "84HZC" = 4, "PKXZL" = 5, "VXY1B" = 6,
  "B9V50" = 1, "ELFHW" = 2, "9ER2B" = 3, "U7DJQ" = 4, "74UOH" = 5, "UWK9R" = 6,
  "LPS" = "", "Neg" = ""
)

Vimentin_summary_df <- Vimentin_summary_df %>%
  mutate(num = mapping[ID])

Vimentin_summary_df <- Vimentin_summary_df %>%
  mutate(Group = fct_relevel(Group, 
                             "PBS Control", 
                             "LPS", 
                             "COVID-19 (-) mTBI (-) ADEs", 
                             "COVID-19 (+) mTBI (-) ADEs", 
                             "COVID-19 (-) mTBI (+) ADEs", 
                             "COVID-19 (+) mTBI (+) ADEs")) %>% arrange(Group, num)

Vimentin_summary_df$ID <- factor(Vimentin_summary_df$ID, levels = Vimentin_summary_df$ID)
Vimentin_summary_df$axis_label <- paste(as.character(Vimentin_summary_df$ID), as.character(Vimentin_summary_df$num), sep="-")
Vimentin_summary_df$axis_label <- factor(Vimentin_summary_df$axis_label, levels = Vimentin_summary_df$axis_label)

# Create a bar graph with facets for group and individual bars for each ID
RA_Ind_Vimentin_plot <- ggplot(Vimentin_summary_df, aes(x = axis_label, y = Mean_vim, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha=0.95) +
  geom_errorbar(aes(ymin = Mean_vim - SE_vim, ymax = Mean_vim + SE_vim), position = position_dodge(width = 0.7), width = 0.25, color="black") +
  theme_minimal()+
  scale_fill_brewer(palette = "Greys") +
  theme(legend.title=element_blank(),legend.position="top",
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black"),
        legend.text = element_text(size = 18, lineheight = .9, family = "sans", face = "bold", colour = "black"))+
  labs(x = "Participant ID",
       y = "Normalized Vimentin MFI",
       fill = "Condition")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Vimentin_summary_df, file.path(output_dir, "figure_4d_RA_Ind_Vimentin_summary.csv"), row.names = FALSE)
write.csv(Vimentin_result_df, file.path(output_dir, "figure_4d_RA_Ind_Vimentin_result.csv"), row.names = FALSE)
ggsave(file.path(tif_dir, "figure_4d_RA_Ind_Vimentin.pdf"), plot = RA_Ind_Vimentin_plot, width = 250, height = 175, units = "mm", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_4d_RA_Ind_Vimentin.pdf")), plot = RA_Ind_Vimentin_plot, width = 250, height = 175, units = "mm", dpi = 600, device = "tiff")
