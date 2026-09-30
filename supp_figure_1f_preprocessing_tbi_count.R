# =============================================================================
# Supplementary Figure 1f preprocessing – Biomarkers by number of mTBIs
# =============================================================================
#
# Description: Accounts for number of mTBIs, not just symptom severity/
#   frequency. qq_tbi_num is heavily right-skewed (range 1-15), so it is
#   categorized as Single (1 mTBI) vs.
#   Multiple (2+ mTBIs) as the primary analysis. Plasma Abeta42/Abeta40 ratio
#   and anti-Abeta1-42 IgG/IgA are compared across these categories
#   (Wilcoxon rank-sum, FDR-adjusted across panels), both pooled ("Overall")
#   and within each mTBI (+) group separately. A continuous Spearman
#   correlation between raw mTBI count and each biomarker is also computed
#   on the same panels, as a robustness check that categorizing mTBI count
#   did not change the conclusion.
#
#   COVID-only participants (no mTBI) are excluded: mTBI count is undefined
#   for them.
#
# Prerequisites: supp_figure_1/supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv;
#   figure_1 output: figure_1i (N3PA Abeta ratio) in figure_1/figure_1_output_data/.
#
# Output: supp_figure_1_tbi_count.pdf, supp_figure_1_tbi_count_continuous.pdf;
#   data in supp_figure_1_output_data
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
# Normalize working directory to repo root regardless of where this is run from.
if (basename(getwd()) == "supp_figure_1") setwd("..")
output_dir <- "supp_figure_1/supp_figure_1_output_data"
input_dir  <- "supp_figure_1/supp_figure_1_input_data"
tif_dir    <- "supp_figure_1"
fig1_out   <- "figure_1/figure_1_output_data"

path_main <- file.path(input_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv")
path_n3pa <- file.path(fig1_out, "figure_1i_N3PA_Abeta_plot_data.csv")
missing <- c(
  if (!file.exists(path_main)) "COAST_Study_Data_Clean_Age_Groups_add_dates.csv (in supp_figure_1/supp_figure_1_input_data)",
  if (!file.exists(path_n3pa)) "figure_1i_N3PA_Abeta_plot_data.csv (run figure_1/figure_1i.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(path_main, stringsAsFactors = FALSE)
N3PA_SRX_Data_Clean_join <- read.csv(path_n3pa, stringsAsFactors = FALSE)

# mTBI (+) groups only; COVID-only (no mTBI) excluded (see header note).
tbi_groups <- c("COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)")

Corr_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::filter(qq_group %in% tbi_groups) %>%
  dplyr::select(participant_id, qq_group, qq_tbi_num, vib_anti_abeta_1_42_igg_iga) %>%
  left_join(N3PA_SRX_Data_Clean_join %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R),
            by = c("participant_id", "qq_group")) %>%
  dplyr::filter(!is.na(qq_tbi_num))

Corr_Data$count_category <- ifelse(Corr_Data$qq_tbi_num >= 2, "Multiple (2+)", "Single (1)")
Corr_Data$count_category <- factor(Corr_Data$count_category, levels = c("Single (1)", "Multiple (2+)"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Corr_Data, file.path(output_dir, "supp_figure_1_tbi_count_input_data.csv"), row.names = FALSE)

cat("Category counts:\n"); print(table(Corr_Data$count_category, Corr_Data$qq_group))

# Long format: one row per participant x biomarker, plus a pooled "Overall" pseudo-group
biomarker_long <- Corr_Data %>%
  tidyr::pivot_longer(cols = c(Abeta42_Abeta40_Ratio_R, vib_anti_abeta_1_42_igg_iga),
                       names_to = "biomarker_name", values_to = "biomarker_value") %>%
  dplyr::mutate(biomarker_name = dplyr::recode(biomarker_name,
                                                Abeta42_Abeta40_Ratio_R = "Aβ42/Aβ40",
                                                vib_anti_abeta_1_42_igg_iga = "Anti-Aβ IgG/IgA")) %>%
  dplyr::filter(!is.na(biomarker_value), !is.na(count_category))

plot_data <- dplyr::bind_rows(
  biomarker_long %>% dplyr::mutate(group_label = "Overall (mTBI +)"),
  biomarker_long %>% dplyr::mutate(group_label = qq_group)
)
plot_data$group_label <- factor(plot_data$group_label,
                                 levels = c("Overall (mTBI +)", tbi_groups))

# Wilcoxon rank-sum per group_label x biomarker panel
compute_panel_stats <- function(df) {
  panels <- df %>% dplyr::distinct(group_label, biomarker_name)
  results <- lapply(seq_len(nrow(panels)), function(i) {
    g <- panels$group_label[i]; b <- panels$biomarker_name[i]
    sub <- df %>% dplyr::filter(group_label == g, biomarker_name == b)
    n_per_cat <- table(droplevels(sub$count_category))
    if (length(n_per_cat) < 2 || any(n_per_cat < 1)) {
      p_value <- NA_real_
    } else {
      p_value <- suppressWarnings(wilcox.test(biomarker_value ~ droplevels(count_category), data = sub)$p.value)
    }
    data.frame(group_label = g, biomarker_name = b, n = nrow(sub),
               n_single = unname(n_per_cat["Single (1)"]), n_multiple = unname(n_per_cat["Multiple (2+)"]),
               wilcoxon_p = p_value)
  })
  do.call(rbind, results)
}
panel_stats <- compute_panel_stats(plot_data)
panel_stats$p_fdr <- p.adjust(panel_stats$wilcoxon_p, method = "fdr")
write.csv(panel_stats, file.path(output_dir, "supp_figure_1_tbi_count_stats.csv"), row.names = FALSE)

Tbi_Count_Plot <- ggplot(plot_data, aes(x = count_category, y = biomarker_value)) +
  geom_boxplot(width = 0.5, fill = "white", color = "black", alpha = 0.7, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.15, height = 0, shape = 21, size = 2.5, stroke = 0.4, fill = "grey60") +
  facet_grid(biomarker_name ~ group_label, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(size = 10, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 10, family = "Arial", colour = "black"),
    axis.title.y = element_text(size = 13, family = "Arial", face = "bold", colour = "black"),
    plot.title = element_text(size = 16, family = "Arial", face = "bold")
  ) +
  labs(
    title = "Biomarkers by Number of mTBIs",
    x = "Number of mTBIs",
    y = "Biomarker Level"
  )

ggsave(file.path(tif_dir, "supp_figure_1_tbi_count.pdf"), Tbi_Count_Plot, width = 12, height = 8, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1_tbi_count.pdf")), Tbi_Count_Plot, width = 12, height = 8, dpi = 600, device = "tiff")

# -----------------------------------------------------------------------
# Continuous companion: Spearman correlation between raw mTBI count and each
# biomarker (same panels as above), as a robustness check that categorizing
# mTBI count did not change the conclusion.
# -----------------------------------------------------------------------
compute_cont_stats <- function(df) {
  panels <- df %>% dplyr::distinct(group_label, biomarker_name)
  results <- lapply(seq_len(nrow(panels)), function(i) {
    g <- panels$group_label[i]; b <- panels$biomarker_name[i]
    sub <- df %>% dplyr::filter(group_label == g, biomarker_name == b)
    ok <- stats::complete.cases(sub$qq_tbi_num, sub$biomarker_value)
    test <- suppressWarnings(stats::cor.test(sub$qq_tbi_num[ok], sub$biomarker_value[ok], method = "spearman"))
    data.frame(group_label = g, biomarker_name = b, n = sum(ok),
               rho = unname(test$estimate), p_value = test$p.value)
  })
  do.call(rbind, results)
}
cont_stats <- compute_cont_stats(plot_data)
cont_stats$p_fdr <- p.adjust(cont_stats$p_value, method = "fdr")
write.csv(cont_stats, file.path(output_dir, "supp_figure_1_tbi_count_continuous_stats.csv"), row.names = FALSE)
cat("\nContinuous (Spearman) companion stats:\n"); print(cont_stats)

Tbi_Count_Continuous_Plot <- ggplot(plot_data, aes(x = qq_tbi_num, y = biomarker_value)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.4, fill = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  facet_grid(biomarker_name ~ group_label, scales = "free") +
  stat_cor(method = "spearman", size = 3.5) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(size = 10, family = "Arial", colour = "black"),
    axis.text.y = element_text(size = 10, family = "Arial", colour = "black"),
    axis.title.y = element_text(size = 13, family = "Arial", face = "bold", colour = "black"),
    plot.title = element_text(size = 16, family = "Arial", face = "bold")
  ) +
  labs(
    title = "Biomarkers vs. Continuous Number of mTBIs",
    x = "Number of mTBIs",
    y = "Biomarker Level"
  )

ggsave(file.path(tif_dir, "supp_figure_1_tbi_count_continuous.pdf"), Tbi_Count_Continuous_Plot, width = 12, height = 8, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1_tbi_count_continuous.pdf")), Tbi_Count_Continuous_Plot, width = 12, height = 8, dpi = 600, device = "tiff")
