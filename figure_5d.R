# =============================================================================
# Figure 5d – HA TLR Inhibitor Vimentin MFI by detailed condition (horizontal bar)
# =============================================================================
#
# Description: Horizontal bar plot of normalized Vimentin MFI by TLR
#   inhibitor/agonist condition (HA astrocyte-derived exosomes). Holm-adjusted
#   t-tests for selected comparisons.
#
# Prerequisites: figure_5_input_data/HA_TLR_Inh_*.csv
#
# Output: figure_5d.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
base_dir   <- if (dir.exists("figure_5")) "figure_5" else "."
input_dir  <- file.path(base_dir, "figure_5_input_data")
output_dir <- file.path(base_dir, "figure_5_output_data")
tif_dir    <- base_dir

library(dplyr)
library(ggplot2)
library(ggsignif)
library(rstatix)
library(tidyr)
library(purrr)

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

# T-tests for selected comparisons (COV_TBI and agonist/inhibitor)
perform_cov_tbi_tests <- function(data) {
  cov_tbi_groups <- data %>% filter(grepl("_COV_TBI$", Condition)) %>% pull(Condition) %>% unique()
  comparisons <- combn(cov_tbi_groups, 2, simplify = FALSE)
  data %>% rstatix::t_test(Normed_Vimentin_MFI ~ Condition, paired = FALSE, comparisons = comparisons) %>%
    adjust_pvalue(method = "holm") %>% add_significance()
}
perform_specific_tests <- function(data, specific_groups) {
  comparisons <- lapply(specific_groups, function(sg) {
    sg_counterpart <- data %>% pull(Condition) %>% unique() %>% grep(paste0("_", sg, "$"), ., value = TRUE)
    if (length(sg_counterpart) == 0) return(NULL)
    list(sg, sg_counterpart)
  })
  comparisons <- comparisons[!sapply(comparisons, is.null)]
  data %>% rstatix::t_test(Normed_Vimentin_MFI ~ Condition, paired = FALSE, comparisons = comparisons) %>%
    adjust_pvalue(method = "holm") %>% add_significance()
}
Vimentin_tests <- perform_cov_tbi_tests(Vimentin_result_df_clean)
specific_groups <- c("Zymosan", "Poly_AU", "LPS", "FLA_BS", "ODN_DSL03")
specific_results <- perform_specific_tests(Vimentin_result_df_clean, specific_groups)
all_results <- bind_rows(Vimentin_tests %>% mutate(test_type = "COV_TBI"), specific_results %>% mutate(test_type = "Specific"))

Vimentin_summary_df_2 <- Vimentin_summary_df_1 %>% filter(Condition != "PEG") %>% filter(!grepl("_Veh$", Condition))
conditions <- c("Exo_COV_TBI", "TL2_C29_COV_TBI", "Zymosan", "TL2_C29_Zymosan", "CU_CPT4a_COV_TBI", "Poly_AU", "CU_CPT4a_Poly_AU", "CLI_095_COV_TBI", "LPS", "CLI_095_LPS", "TH1020_COV_TBI", "FLA_BS", "TH1020_FLA_BS", "ODN_INH18_COV_TBI", "ODN_DSL03", "ODN_INH18_ODN_DSL03")
Vimentin_summary_df_2 <- Vimentin_summary_df_2 %>% filter(Condition %in% conditions)

replace_substring <- function(x) {
  x <- gsub("_COV_TBI$", "-COVID-19 (+) mTBI (+) ADEs", x)
  x <- gsub("^Exo-", "", x)
  x <- gsub("_", "-", x)
  x
}
Vimentin_summary_df_2$Condition <- sapply(Vimentin_summary_df_2$Condition, replace_substring)
Vimentin_summary_df_2$Condition <- factor(Vimentin_summary_df_2$Condition,
  levels = c("COVID-19 (+) mTBI (+) ADEs", "ODN-DSL03", "ODN-INH18-ODN-DSL03", "ODN-INH18-COVID-19 (+) mTBI (+) ADEs", "FLA-BS", "TH1020-FLA-BS", "TH1020-COVID-19 (+) mTBI (+) ADEs", "LPS", "CLI-095-LPS", "CLI-095-COVID-19 (+) mTBI (+) ADEs", "Poly-AU", "CU-CPT4a-Poly-AU", "CU-CPT4a-COVID-19 (+) mTBI (+) ADEs", "Zymosan", "TL2-C29-Zymosan", "TL2-C29-COVID-19 (+) mTBI (+) ADEs"))
Vimentin_summary_df_2_clean <- Vimentin_summary_df_2 %>% filter(Condition != "NA")

format_pvalue_4digits <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
get_formatted_pvalue <- function(index) format_pvalue_4digits(all_results$p.adj[index])

max_y <- max(Vimentin_summary_df_2_clean$Mean + Vimentin_summary_df_2_clean$SE, na.rm = TRUE)

p <- ggplot(Vimentin_summary_df_2_clean, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.95) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.7), width = 0.25, color = "black") +
  theme_minimal() + scale_fill_grey(start = 0.8, end = 0.2) +
  theme(legend.title = element_blank(), legend.position = "none", plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid.major = element_line(color = "grey90"), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
        axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
        legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black")) +
  labs(x = "Condition", y = "Normalized Vimentin MFI") + coord_flip() + scale_y_continuous(limits = c(0, max_y * 2.3))

p <- p +
  geom_signif(comparisons = list(c("ODN-DSL03", "ODN-INH18-ODN-DSL03")), annotations = get_formatted_pvalue(20), y_position = max_y * 0.9, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("FLA-BS", "TH1020-FLA-BS")), annotations = get_formatted_pvalue(19), y_position = max_y * 1.0, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("LPS", "CLI-095-LPS")), annotations = get_formatted_pvalue(18), y_position = max_y * 1.1, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("Poly-AU", "CU-CPT4a-Poly-AU")), annotations = get_formatted_pvalue(17), y_position = max_y * 1.2, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("Zymosan", "TL2-C29-Zymosan")), annotations = get_formatted_pvalue(16), y_position = max_y * 1.3, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("CLI-095-COVID-19 (+) mTBI (+) ADEs", "CU-CPT4a-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(1), y_position = max_y * 1.4, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("CLI-095-COVID-19 (+) mTBI (+) ADEs", "ODN-INH18-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(3), y_position = max_y * 1.5, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("CU-CPT4a-COVID-19 (+) mTBI (+) ADEs", "COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(6), y_position = max_y * 1.6, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("CU-CPT4a-COVID-19 (+) mTBI (+) ADEs", "ODN-INH18-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(7), y_position = max_y * 1.7, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("CU-CPT4a-COVID-19 (+) mTBI (+) ADEs", "TH1020-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(8), y_position = max_y * 1.8, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("CU-CPT4a-COVID-19 (+) mTBI (+) ADEs", "TL2-C29-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(9), y_position = max_y * 1.9, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("COVID-19 (+) mTBI (+) ADEs", "ODN-INH18-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(10), y_position = max_y * 2, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("ODN-INH18-COVID-19 (+) mTBI (+) ADEs", "TH1020-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(13), y_position = max_y * 2.1, tip_length = 0, size = 1, textsize = 8) +
  geom_signif(comparisons = list(c("ODN-INH18-COVID-19 (+) mTBI (+) ADEs", "TL2-C29-COVID-19 (+) mTBI (+) ADEs")), annotations = get_formatted_pvalue(14), y_position = max_y * 2.2, tip_length = 0, size = 1, textsize = 8)

print(p)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Vimentin_summary_df_2_clean, file.path(output_dir, "figure_5d_Vimentin_HA_TLR_Inh_summary.csv"), row.names = FALSE)
write.csv(Vimentin_result_df_clean, file.path(output_dir, "figure_5d_Vimentin_HA_TLR_Inh_result.csv"), row.names = FALSE)
ggsave(file.path(tif_dir, "figure_5d.pdf"), plot = p, width = 400, height = 200, units = "mm", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_5d.pdf")), plot = p, width = 400, height = 200, units = "mm", dpi = 600, device = "tiff")
