# =============================================================================
# Supplementary Figure 1h preprocessing – COVID-attributed outcomes by number of prior mTBIs
# =============================================================================
#
# Description: Companion to supp_figure_1_preprocessing_mtbi_covid_gap.R,
#   which tests whether the timing of a participant's most recent prior
#   mTBI relates to COVID-19 outcomes. This script instead checks how many
#   prior mTBIs someone had (Single (1) vs. Multiple (2+), the same
#   categorization supp_figure_1f_preprocessing_tbi_count.R uses for
#   overall mTBI burden), restricted to the double-exposed group.
#
#   Restricted to participants with exactly one COVID-19 infection: the
#   COVID-19 severity/frequency scores sum a change score across every
#   COVID-19 infection a participant tested positive for, so including
#   participants with multiple infections would make the outcome ambiguous.
#
# Prerequisites: supp_figure_1/supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv;
#   supp_figure_1/supp_figure_1_input_data/mtbi_covid_gap.csv (run
#   supp_figure_1_preprocessing_mtbi_covid_gap.R first);
#   figure_1 output: figure_1i (N3PA Abeta ratio) in figure_1/figure_1_output_data/
#   (or supp_figure_1/supp_figure_1_input_data/figure_1i_N3PA_Abeta_plot_data.csv).
#
# Output: supp_figure_1_mtbi_count.pdf, supp_figure_1_mtbi_count_continuous.pdf;
#   data in supp_figure_1_output_data
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
if (basename(getwd()) == "supp_figure_1") setwd("..")
output_dir <- "supp_figure_1/supp_figure_1_output_data"
input_dir  <- "supp_figure_1/supp_figure_1_input_data"
tif_dir    <- "supp_figure_1"
fig1_out   <- "figure_1/figure_1_output_data"

path_main <- file.path(input_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv")
path_gap  <- file.path(input_dir, "mtbi_covid_gap.csv")
path_n3pa <- file.path(fig1_out, "figure_1i_N3PA_Abeta_plot_data.csv")
if (!file.exists(path_n3pa)) path_n3pa <- file.path(input_dir, "figure_1i_N3PA_Abeta_plot_data.csv")
missing <- c(
  if (!file.exists(path_main)) "COAST_Study_Data_Clean_Age_Groups_add_dates.csv (in supp_figure_1/supp_figure_1_input_data)",
  if (!file.exists(path_gap))  "mtbi_covid_gap.csv (run supp_figure_1_preprocessing_mtbi_covid_gap.R first)",
  if (!file.exists(path_n3pa)) "figure_1i_N3PA_Abeta_plot_data.csv (run figure_1/figure_1i.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(path_main, stringsAsFactors = FALSE)
Gap_Data <- read.csv(path_gap, stringsAsFactors = FALSE)
N3PA_SRX_Data_Clean_join <- read.csv(path_n3pa, stringsAsFactors = FALSE)

# Same double-exposed population as supp_figure_1_mtbi_covid_gap.R: participants
# with at least one mTBI before their first COVID-19 infection (n_tbi_before_covid
# defined). n_tbi_before_covid is heavily right-skewed (1-10), so it is
# categorized as Single (1) vs. Multiple (2+), matching the panel f script's
# supp_figure_1_tbi_count.R convention.
Corr_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::filter(qq_group == "COVID-19 (+) mTBI (+)") %>%
  dplyr::select(participant_id, qq_group, covid_pos_test_num, recode_covid_neuro_psych_sev_score,
                recode_covid_neuro_psych_freq_score, vib_anti_abeta_1_42_igg_iga) %>%
  left_join(Gap_Data %>% dplyr::select(participant_id, n_tbi_before_covid), by = "participant_id") %>%
  left_join(N3PA_SRX_Data_Clean_join %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R),
            by = c("participant_id", "qq_group")) %>%
  # n_tbi_before_covid == 0 means this participant's only mTBI(s) are on/after
  # COVID-19 (no prior mTBI at all) - excluded here, same as the gap analysis,
  # since "Single vs. Multiple prior mTBI" presupposes at least one.
  dplyr::filter(!is.na(n_tbi_before_covid), n_tbi_before_covid >= 1) %>%
  # Single COVID-19 infection only (see header note): recode_covid_neuro_psych_*
  # sums across all infections, so participants with 2+ would blend outcomes
  # from an infection unrelated to the mTBI count being tested here.
  dplyr::filter(covid_pos_test_num == 1)

cat("Double-positive participants with >=1 pre-COVID-19 mTBI and a single COVID-19 infection: n =", nrow(Corr_Data), "\n")

Corr_Data$count_category <- ifelse(Corr_Data$n_tbi_before_covid >= 2, "Multiple (2+)", "Single (1)")
Corr_Data$count_category <- factor(Corr_Data$count_category, levels = c("Single (1)", "Multiple (2+)"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Corr_Data, file.path(output_dir, "supp_figure_1_mtbi_count_input_data.csv"), row.names = FALSE)

cat("\nCategory counts:\n"); print(table(Corr_Data$count_category))

outcome_long <- Corr_Data %>%
  tidyr::pivot_longer(cols = c(recode_covid_neuro_psych_sev_score, recode_covid_neuro_psych_freq_score,
                                Abeta42_Abeta40_Ratio_R, vib_anti_abeta_1_42_igg_iga),
                       names_to = "outcome_name", values_to = "outcome_value") %>%
  dplyr::mutate(outcome_name = dplyr::recode(outcome_name,
                                              recode_covid_neuro_psych_sev_score = "COVID Severity",
                                              recode_covid_neuro_psych_freq_score = "COVID Frequency",
                                              Abeta42_Abeta40_Ratio_R = "Aβ42/Aβ40",
                                              vib_anti_abeta_1_42_igg_iga = "Anti-Aβ IgG/IgA")) %>%
  dplyr::filter(!is.na(outcome_value), !is.na(count_category))
outcome_long$outcome_name <- factor(outcome_long$outcome_name,
                                     levels = c("COVID Severity", "COVID Frequency", "Aβ42/Aβ40", "Anti-Aβ IgG/IgA"))

# Wilcoxon rank-sum (2 groups) per outcome, FDR across outcomes.
compute_panel_stats <- function(df) {
  outcomes <- unique(df$outcome_name)
  results <- lapply(outcomes, function(o) {
    sub <- df %>% dplyr::filter(outcome_name == o)
    n_per_cat <- table(droplevels(sub$count_category))
    p_value <- if (length(n_per_cat) < 2 || any(n_per_cat < 1)) {
      NA_real_
    } else {
      suppressWarnings(wilcox.test(outcome_value ~ droplevels(count_category), data = sub)$p.value)
    }
    data.frame(outcome_name = as.character(o), n = nrow(sub),
               n_single = unname(n_per_cat["Single (1)"]), n_multiple = unname(n_per_cat["Multiple (2+)"]),
               wilcoxon_p = p_value)
  })
  do.call(rbind, results)
}
panel_stats <- compute_panel_stats(outcome_long)
panel_stats$p_fdr <- p.adjust(panel_stats$wilcoxon_p, method = "fdr")
write.csv(panel_stats, file.path(output_dir, "supp_figure_1_mtbi_count_stats.csv"), row.names = FALSE)
cat("\nWilcoxon (Single vs. Multiple prior mTBI) stats:\n"); print(panel_stats, row.names = FALSE)

# Same per-outcome labeled y-axis / explicit-margin label placement as
# supp_figure_1_mtbi_covid_gap.R (see that script for why: ggpubr's "top" npc
# position is computed from the data range, not the expanded axis range, so
# it can still land on an extreme point).
outcome_meta <- data.frame(
  outcome_name = c("COVID Severity", "COVID Frequency", "Aβ42/Aβ40", "Anti-Aβ IgG/IgA"),
  y_label = c("COVID-19 Symptom\nSeverity Score", "COVID-19 Symptom\nFrequency Score",
              "Plasma Aβ42/Aβ40 Ratio", "Anti-Aβ1-42 (IgG/IgA)"),
  stringsAsFactors = FALSE
)

build_count_panel <- function(outcome) {
  df <- outcome_long %>% dplyr::filter(outcome_name == outcome)
  y_lab <- outcome_meta$y_label[outcome_meta$outcome_name == outcome]
  y_range <- range(df$outcome_value, na.rm = TRUE)
  label_y <- y_range[2] + 0.12 * diff(y_range)
  ggplot(df, aes(x = count_category, y = outcome_value)) +
    geom_boxplot(width = 0.5, fill = "white", color = "black", alpha = 0.7, outlier.shape = NA, linewidth = 0.6) +
    geom_jitter(width = 0.15, height = 0, shape = 21, size = 2.5, stroke = 0.4, fill = "grey60") +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5, label.y = label_y) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, family = "Arial", face = "bold", colour = "black", hjust = 0.5),
      axis.text.x = element_text(size = 9, family = "Arial", face = "bold", colour = "black"),
      axis.text.y = element_text(size = 10, family = "Arial", colour = "black"),
      axis.title.y = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
      axis.title.x = element_blank()
    ) +
    labs(title = outcome, y = y_lab)
}
count_panels <- lapply(outcome_meta$outcome_name, build_count_panel)

Count_Plot <- ggpubr::ggarrange(plotlist = count_panels, nrow = 1)
Count_Plot <- ggpubr::annotate_figure(
  Count_Plot,
  top = ggpubr::text_grob("COVID-19 Outcomes by Number of Prior mTBIs, Double-Exposed Group, Single COVID-19 Infection",
                           face = "bold", size = 15, family = "Arial"),
  bottom = ggpubr::text_grob("Number of mTBIs Before First COVID-19 Infection",
                              face = "bold", size = 12, family = "Arial")
)
ggsave(file.path(tif_dir, "supp_figure_1_mtbi_count.pdf"), Count_Plot, width = 14, height = 6.5, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1_mtbi_count.pdf")), Count_Plot, width = 14, height = 6.5, dpi = 600, device = "tiff")

# -----------------------------------------------------------------------
# Continuous companion: Spearman correlation between raw pre-COVID-19 mTBI
# count and each outcome, as a robustness check that categorizing count did
# not change the conclusion (same convention as supp_figure_1_mtbi_covid_gap.R).
# -----------------------------------------------------------------------
compute_cont_stats <- function(df) {
  outcomes <- unique(df$outcome_name)
  results <- lapply(outcomes, function(o) {
    sub <- df %>% dplyr::filter(outcome_name == o)
    ok <- stats::complete.cases(sub$n_tbi_before_covid, sub$outcome_value)
    test <- suppressWarnings(stats::cor.test(sub$n_tbi_before_covid[ok], sub$outcome_value[ok], method = "spearman"))
    data.frame(outcome_name = as.character(o), n = sum(ok), rho = unname(test$estimate), p_value = test$p.value)
  })
  do.call(rbind, results)
}
cont_stats <- compute_cont_stats(outcome_long)
cont_stats$p_fdr <- p.adjust(cont_stats$p_value, method = "fdr")
write.csv(cont_stats, file.path(output_dir, "supp_figure_1_mtbi_count_continuous_stats.csv"), row.names = FALSE)
cat("\nContinuous (Spearman) companion stats:\n"); print(cont_stats, row.names = FALSE)

build_count_continuous_panel <- function(outcome) {
  df <- outcome_long %>% dplyr::filter(outcome_name == outcome)
  y_lab <- outcome_meta$y_label[outcome_meta$outcome_name == outcome]
  y_range <- range(df$outcome_value, na.rm = TRUE)
  label_y <- y_range[2] + 0.12 * diff(y_range)
  label_x <- min(df$n_tbi_before_covid, na.rm = TRUE)
  ggplot(df, aes(x = n_tbi_before_covid, y = outcome_value)) +
    geom_point(shape = 21, size = 2.5, stroke = 0.4, fill = "grey60") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    stat_cor(method = "spearman", size = 3.5, label.x = label_x, label.y = label_y) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, family = "Arial", face = "bold", colour = "black", hjust = 0.5),
      axis.text.x = element_text(size = 10, family = "Arial", colour = "black"),
      axis.text.y = element_text(size = 10, family = "Arial", colour = "black"),
      axis.title.y = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
      axis.title.x = element_blank()
    ) +
    labs(title = outcome, y = y_lab)
}
count_continuous_panels <- lapply(outcome_meta$outcome_name, build_count_continuous_panel)

Count_Continuous_Plot <- ggpubr::ggarrange(plotlist = count_continuous_panels, nrow = 1)
Count_Continuous_Plot <- ggpubr::annotate_figure(
  Count_Continuous_Plot,
  top = ggpubr::text_grob("COVID-19 Outcomes vs. Continuous Number of Prior mTBIs, Double-Exposed Group, Single COVID-19 Infection",
                           face = "bold", size = 15, family = "Arial"),
  bottom = ggpubr::text_grob("Number of mTBIs Before First COVID-19 Infection",
                              face = "bold", size = 12, family = "Arial")
)
ggsave(file.path(tif_dir, "supp_figure_1_mtbi_count_continuous.pdf"), Count_Continuous_Plot, width = 14, height = 6.5, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1_mtbi_count_continuous.pdf")), Count_Continuous_Plot, width = 14, height = 6.5, dpi = 600, device = "tiff")
