# =============================================================================
# Supplementary Figure 1i preprocessing – COVID-attributed outcomes by
#   severity/frequency of the specific prior mTBI
# =============================================================================
#
# Description: Tests whether the severity of a participant's specific
#   prior mTBI (not just its timing or count) relates to their COVID-19
#   outcomes, restricted to the double-exposed group. Uses a severity/
#   frequency score restricted to that one target mTBI incidence (from
#   supp_figure_1_preprocessing_mtbi_severity_target.R), since the
#   whole-life aggregate score sums across every mTBI a participant ever
#   had, which would conflate unrelated injuries.
#
#   Only the two COVID-19 symptom outcomes are plotted; plasma Aβ42/Aβ40
#   and anti-Aβ1-42 IgG/IgA were also tested against both predictors (all
#   p > 0.4, see supp_figure_1_output_data/supp_figure_1_mtbi_*_continuous_stats.csv)
#   but are not shown since neither relationship was significant.
#
# Prerequisites: supp_figure_1/supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv;
#   supp_figure_1/supp_figure_1_input_data/mtbi_severity_target_incidence.csv
#   (run supp_figure_1_preprocessing_mtbi_severity_target.R first).
#
# Output: supp_figure_1_mtbi_severity_continuous.pdf, supp_figure_1_mtbi_frequency_continuous.pdf;
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
path_sev  <- file.path(input_dir, "mtbi_severity_target_incidence.csv")
path_n3pa <- file.path(fig1_out, "figure_1i_N3PA_Abeta_plot_data.csv")
if (!file.exists(path_n3pa)) path_n3pa <- file.path(input_dir, "figure_1i_N3PA_Abeta_plot_data.csv")
missing <- c(
  if (!file.exists(path_main)) "COAST_Study_Data_Clean_Age_Groups_add_dates.csv (in supp_figure_1/supp_figure_1_input_data)",
  if (!file.exists(path_sev))  "mtbi_severity_target_incidence.csv (run supp_figure_1_mtbi_severity_preprocessing.R first)",
  if (!file.exists(path_n3pa)) "figure_1i_N3PA_Abeta_plot_data.csv (run figure_1/figure_1i.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dunn.test)
library(scales)

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(path_main, stringsAsFactors = FALSE)
Severity_Data <- read.csv(path_sev, stringsAsFactors = FALSE)
N3PA_SRX_Data_Clean_join <- read.csv(path_n3pa, stringsAsFactors = FALSE)

Corr_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::filter(qq_group == "COVID-19 (+) mTBI (+)") %>%
  dplyr::select(participant_id, qq_group, covid_pos_test_num, recode_covid_neuro_psych_sev_score,
                recode_covid_neuro_psych_freq_score, vib_anti_abeta_1_42_igg_iga) %>%
  left_join(Severity_Data %>% dplyr::select(participant_id, target_incidence,
                                             tbi_severity_target_incidence, tbi_frequency_target_incidence),
            by = "participant_id") %>%
  left_join(N3PA_SRX_Data_Clean_join %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R),
            by = c("participant_id", "qq_group")) %>%
  dplyr::filter(!is.na(target_incidence)) %>%
  # Single COVID-19 infection only (see header note).
  dplyr::filter(covid_pos_test_num == 1)

cat("Double-positive participants with a target incidence and a single COVID-19 infection: n =", nrow(Corr_Data), "\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Corr_Data, file.path(output_dir, "supp_figure_1_mtbi_severity_input_data.csv"), row.names = FALSE)

# All 4 outcomes are still tested and saved to the stats CSVs (full record,
# including the non-significant biomarkers); only the 2 COVID-19 symptom
# outcomes are plotted below.
outcome_meta <- data.frame(
  outcome_name = c("COVID Severity", "COVID Frequency", "Aβ42/Aβ40", "Anti-Aβ IgG/IgA"),
  y_label = c("COVID-19 Symptom\nSeverity Score", "COVID-19 Symptom\nFrequency Score",
              "Plasma Aβ42/Aβ40 Ratio", "Anti-Aβ1-42 (IgG/IgA)"),
  stringsAsFactors = FALSE
)
plotted_outcomes <- c("COVID Severity", "COVID Frequency")

make_outcome_long <- function(df) {
  long <- df %>%
    tidyr::pivot_longer(cols = c(recode_covid_neuro_psych_sev_score, recode_covid_neuro_psych_freq_score,
                                  Abeta42_Abeta40_Ratio_R, vib_anti_abeta_1_42_igg_iga),
                         names_to = "outcome_name", values_to = "outcome_value") %>%
    dplyr::mutate(outcome_name = dplyr::recode(outcome_name,
                                                recode_covid_neuro_psych_sev_score = "COVID Severity",
                                                recode_covid_neuro_psych_freq_score = "COVID Frequency",
                                                Abeta42_Abeta40_Ratio_R = "Aβ42/Aβ40",
                                                vib_anti_abeta_1_42_igg_iga = "Anti-Aβ IgG/IgA")) %>%
    dplyr::filter(!is.na(outcome_value))
  long$outcome_name <- factor(long$outcome_name,
                               levels = c("COVID Severity", "COVID Frequency", "Aβ42/Aβ40", "Anti-Aβ IgG/IgA"))
  long
}

# Kruskal-Wallis (omnibus) + Dunn post-hoc (Holm) per outcome, for a given
# tertile-categorized predictor column.
compute_panel_stats <- function(df, category_col) {
  outcomes <- unique(df$outcome_name)
  results <- lapply(outcomes, function(o) {
    sub <- df %>% dplyr::filter(outcome_name == o, !is.na(.data[[category_col]]))
    kw <- kruskal.test(sub$outcome_value ~ droplevels(sub[[category_col]]))
    dunn <- suppressWarnings(dunn.test(sub$outcome_value, sub[[category_col]], method = "holm", table = FALSE, list = FALSE))
    dunn_df <- data.frame(comparison = dunn$comparisons, Z = dunn$Z, p_value = dunn$P, p_holm = dunn$P.adjusted)
    dunn_df$outcome_name <- as.character(o)
    dunn_df$kruskal_p <- kw$p.value
    dunn_df$n <- nrow(sub)
    dunn_df
  })
  do.call(rbind, results)
}

compute_cont_stats <- function(df, predictor_col) {
  outcomes <- unique(df$outcome_name)
  results <- lapply(outcomes, function(o) {
    sub <- df %>% dplyr::filter(outcome_name == o)
    ok <- stats::complete.cases(sub[[predictor_col]], sub$outcome_value)
    test <- suppressWarnings(stats::cor.test(sub[[predictor_col]][ok], sub$outcome_value[ok], method = "spearman"))
    data.frame(outcome_name = as.character(o), n = sum(ok), rho = unname(test$estimate), p_value = test$p.value)
  })
  do.call(rbind, results)
}

# Styling matches this manuscript's own figure_1/efigure_1 convention (see
# figure_1i.R): theme_minimal with the gridlines stripped out entirely in
# favor of solid black axis lines/ticks, bold Arial text, and the same
# burnt-orange accent (#DD4726) already used there for the highlighted
# element - here, the trend line. Point color is a dedicated steel blue
# (distinct from the accent and from the manuscript's own grayscale point
# fills) so the two-color pairing (points vs. trend) reads cleanly at a
# glance.
point_fill   <- "#2E6E9E"
accent_color <- "#DD4726"

manuscript_theme <- theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(size = 15, family = "Arial", face = "bold", colour = "black", hjust = 0.5),
    axis.text.x = element_text(size = 13, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 13, family = "Arial", face = "bold", colour = "black"),
    axis.title.y = element_text(size = 14, family = "Arial", face = "bold", colour = "black", margin = margin(r = 12)),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 10, r = 16, b = 6, l = 6)
  )

# Prior mTBI severity/frequency are heavily right-skewed with a genuine floor
# at 0 (unlike the gap variable, ~25% of participants score exactly 0), so a
# true log scale isn't usable. scales::pseudo_log_trans() behaves linearly
# near 0 and log-like further out, spreading the clustered majority of the
# sample out without dropping the zeros - Spearman's rho/p (rank-based) are
# unaffected by the transform.
build_continuous_panel <- function(outcome_long, predictor_col, outcome, x_breaks) {
  df <- outcome_long %>% dplyr::filter(outcome_name == outcome, !is.na(.data[[predictor_col]]))
  y_lab <- outcome_meta$y_label[outcome_meta$outcome_name == outcome]
  y_range <- range(df$outcome_value, na.rm = TRUE)
  label_y <- y_range[2] + 0.12 * diff(y_range)
  label_x <- min(df[[predictor_col]], na.rm = TRUE)
  ggplot(df, aes(x = .data[[predictor_col]], y = outcome_value)) +
    geom_point(shape = 21, size = 4, stroke = 0.6, colour = "black", fill = point_fill, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = accent_color, fill = accent_color, alpha = 0.15, linewidth = 1) +
    stat_cor(method = "spearman", size = 4.2, label.x = label_x, label.y = label_y, fontface = "bold") +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = x_breaks) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    manuscript_theme +
    labs(title = outcome, y = y_lab)
}

# Runs the analysis for one predictor (tbi_severity_target_incidence or
# tbi_frequency_target_incidence). All 4 outcomes are tested and saved to the
# stats CSVs (tertile Kruskal-Wallis/Dunn AND continuous Spearman), but only
# COVID Severity/Frequency are plotted:
#   - The tertile (Low/Moderate/High) categorical comparison forces this
#     heavily right-skewed variable into equal-COUNT tertiles, producing a
#     "High" bucket spanning nearly an order of magnitude (e.g. 11-85 for
#     frequency) that dilutes exactly the signal the continuous view shows
#     cleanly - so it's kept as a supplementary stat, not plotted.
#   - Aβ42/Aβ40 and Anti-Aβ IgG/IgA are not significant against either
#     predictor (all p > 0.4; see the continuous stats CSV) and are dropped
#     from the figure.
run_predictor_analysis <- function(predictor_col, predictor_label, file_stub, title_word, x_breaks) {
  df <- Corr_Data %>% dplyr::filter(!is.na(.data[[predictor_col]]))
  cuts <- quantile(df[[predictor_col]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  category_col <- paste0(predictor_col, "_category")
  df[[category_col]] <- cut(df[[predictor_col]], breaks = cuts, include.lowest = TRUE,
                             labels = c("Low", "Moderate", "High"))

  cat("\n---", title_word, "---\n")
  cat("Tertile cutpoints (supplementary stats only, not plotted):\n"); print(cuts)
  cat("Category counts:\n"); print(table(df[[category_col]]))

  outcome_long <- make_outcome_long(df)

  panel_stats <- compute_panel_stats(outcome_long, category_col)
  kruskal_p_table <- unique(panel_stats[, c("outcome_name", "kruskal_p")])
  kruskal_p_table$kruskal_p_fdr <- p.adjust(kruskal_p_table$kruskal_p, method = "fdr")
  panel_stats <- merge(panel_stats, kruskal_p_table, by = c("outcome_name", "kruskal_p"))
  panel_stats <- panel_stats[, c("outcome_name", "n", "kruskal_p", "kruskal_p_fdr", "comparison", "Z", "p_value", "p_holm")]
  write.csv(panel_stats, file.path(output_dir, paste0(file_stub, "_stats.csv")), row.names = FALSE)
  cat("Tertile (Kruskal-Wallis / Dunn) stats (all 4 outcomes):\n"); print(panel_stats, row.names = FALSE)

  cont_stats <- compute_cont_stats(outcome_long, predictor_col)
  cont_stats$p_fdr <- p.adjust(cont_stats$p_value, method = "fdr")
  write.csv(cont_stats, file.path(output_dir, paste0(file_stub, "_continuous_stats.csv")), row.names = FALSE)
  cat("\nContinuous (Spearman) companion stats (all 4 outcomes):\n"); print(cont_stats, row.names = FALSE)

  continuous_panels <- lapply(plotted_outcomes, function(o) build_continuous_panel(outcome_long, predictor_col, o, x_breaks))
  continuous_plot <- ggpubr::ggarrange(plotlist = continuous_panels, nrow = 1)
  continuous_plot <- ggpubr::annotate_figure(
    continuous_plot,
    bottom = ggpubr::text_grob(predictor_label, face = "bold", size = 14, family = "Arial")
  )
  ggsave(file.path(tif_dir, paste0(file_stub, "_continuous.pdf")), continuous_plot, width = 10, height = 6, dpi = 600, device = grDevices::cairo_pdf)
  ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, paste0(file_stub, "_continuous.pdf"))), continuous_plot, width = 10, height = 6, dpi = 600, device = "tiff")
}

run_predictor_analysis("tbi_severity_target_incidence",
                        "Severity of the Most Recent Prior mTBI (Before COVID-19)",
                        "supp_figure_1_mtbi_severity", "Severity", x_breaks = c(0, 2, 5, 10, 20, 40, 75))
run_predictor_analysis("tbi_frequency_target_incidence",
                        "Frequency of Symptoms from the Most Recent Prior mTBI (Before COVID-19)",
                        "supp_figure_1_mtbi_frequency", "Frequency", x_breaks = c(0, 2, 5, 10, 20, 40, 85))
