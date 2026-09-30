# =============================================================================
# Supplementary Figure 1i – Does mTBI symptom severity/frequency predict
#   COVID-19 outcomes after controlling for other variables?
# =============================================================================
#
# Description: supp_figure_1i_preprocessing_mtbi_severity.R reports
# unadjusted (bivariate) Spearman correlations between the severity/
# frequency of a participant's target prior mTBI and their COVID-19
# symptom severity/frequency. This script adjusts for n_tbi_before_covid
# (number of mTBIs before the first COVID-19 infection, computed by
# supp_figure_1_preprocessing_mtbi_covid_gap.R) as "# of past injuries."
# # of past COVID-19 infections is not included: this analysis is already
# restricted to participants with exactly one, so it has zero variance in
# this sample and cannot be used as a covariate.
#
#   n = 40 supports at most ~2-3 predictors before the model becomes
#   unstable (rule of thumb ~10-15 observations per parameter), so the
#   model is kept to one main predictor (mTBI severity OR frequency) plus
#   this covariate. Both continuous variables are rank-transformed before
#   fitting (rank(outcome) ~ rank(predictor) + rank(injury count)),
#   preserving the same robustness to skew/outliers as the Spearman
#   correlations reported elsewhere in this folder.
#
# Prerequisites: supp_figure_1/supp_figure_1_output_data/supp_figure_1_mtbi_severity_input_data.csv
#   (run supp_figure_1i_preprocessing_mtbi_severity.R first);
#   supp_figure_1/supp_figure_1_input_data/mtbi_covid_gap.csv (for n_tbi_before_covid).
#
# Output: supp_figure_1i_mtbi_severity_regression_forest.pdf; adjusted model
#   coefficients in supp_figure_1_output_data/supp_figure_1i_mtbi_severity_regression_stats.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) == "supp_figure_1") setwd("..")
output_dir <- "supp_figure_1/supp_figure_1_output_data"
input_dir  <- "supp_figure_1/supp_figure_1_input_data"
tif_dir    <- "supp_figure_1"
fig1_out   <- "figure_1/figure_1_output_data"

path_input <- file.path(output_dir, "supp_figure_1_mtbi_severity_input_data.csv")
path_gap   <- file.path(input_dir, "mtbi_covid_gap.csv")
missing <- c(
  if (!file.exists(path_input)) "supp_figure_1_mtbi_severity_input_data.csv (run supp_figure_1i_preprocessing_mtbi_severity.R first)",
  if (!file.exists(path_gap))   "mtbi_covid_gap.csv (run supp_figure_1_preprocessing_mtbi_covid_gap.R first)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(ggplot2)
library(car)  # vif()
library(ggh4x)  # facet_grid2() for per-panel y-axis dropping

Corr_Data <- read.csv(path_input, stringsAsFactors = FALSE)
gap_data  <- read.csv(path_gap, stringsAsFactors = FALSE)

Corr_Data <- Corr_Data %>%
  dplyr::left_join(gap_data %>% dplyr::select(participant_id, n_tbi_before_covid), by = "participant_id")

cat("n =", nrow(Corr_Data), "; complete injury count:",
    sum(!is.na(Corr_Data$n_tbi_before_covid)), "\n")

# Multicollinearity check among candidate predictors (Spearman, since these
# are the same skewed variables used throughout this folder).
cat("\nSpearman correlations among candidate predictors (checking for collinearity before adjusting):\n")
print(round(cor(Corr_Data[, c("tbi_severity_target_incidence", "tbi_frequency_target_incidence",
                               "n_tbi_before_covid")],
                 method = "spearman", use = "pairwise.complete.obs"), 2))

predictor_specs <- list(
  list(col = "tbi_severity_target_incidence",  label = "Severity"),
  list(col = "tbi_frequency_target_incidence", label = "Frequency")
)
outcome_specs <- list(
  list(col = "recode_covid_neuro_psych_sev_score",  label = "COVID Severity"),
  list(col = "recode_covid_neuro_psych_freq_score", label = "COVID Frequency")
)

# Fits rank(outcome) ~ rank(predictor) + rank(injury count) and returns BOTH
# terms - the predictor of interest AND the adjustment covariate - so the
# full model is visible, not just the term we're testing.
fit_adjusted_model <- function(df, predictor_col, predictor_name, outcome_col) {
  sub <- df %>%
    dplyr::select(all_of(c(predictor_col, outcome_col, "n_tbi_before_covid"))) %>%
    stats::na.omit()
  sub_ranked <- as.data.frame(lapply(sub, rank))
  form <- stats::as.formula(paste0(outcome_col, " ~ ", predictor_col, " + n_tbi_before_covid"))
  model <- stats::lm(form, data = sub_ranked)
  vifs <- car::vif(model)
  coefs <- summary(model)$coefficients
  ci <- suppressMessages(stats::confint(model))
  term_row <- function(term, term_name, role) {
    data.frame(
      term_name = term_name, term_role = role, n = nrow(sub),
      estimate = coefs[term, "Estimate"], ci_lower = ci[term, 1], ci_upper = ci[term, 2],
      p_value = coefs[term, "Pr(>|t|)"], r_squared = summary(model)$r.squared, vif = unname(vifs[term])
    )
  }
  rbind(
    term_row(predictor_col, predictor_name, "Predictor of interest"),
    term_row("n_tbi_before_covid", "Number of prior mTBIs", "Adjustment covariate")
  )
}

results <- do.call(rbind, lapply(predictor_specs, function(p) {
  do.call(rbind, lapply(outcome_specs, function(o) {
    res <- fit_adjusted_model(Corr_Data, p$col, paste0("mTBI symptom ", tolower(p$label)), o$col)
    res$predictor_label <- p$label
    res$outcome_label <- o$label
    res
  }))
}))
# FDR correction applies only to the predictor-of-interest tests (the actual
# hypotheses being evaluated, matching the original 4-test scope) - the
# covariate's own p-value is shown for transparency, not as a formal test.
results$p_fdr <- NA_real_
is_predictor <- results$term_role == "Predictor of interest"
results$p_fdr[is_predictor] <- p.adjust(results$p_value[is_predictor], method = "fdr")
results <- results[, c("predictor_label", "outcome_label", "term_role", "term_name", "n",
                        "estimate", "ci_lower", "ci_upper", "p_value", "p_fdr", "r_squared", "vif")]

write.csv(results, file.path(output_dir, "supp_figure_1i_mtbi_severity_regression_stats.csv"), row.names = FALSE)
cat("\nAdjusted (rank-regression) effect of mTBI severity/frequency on COVID outcomes,\n",
    "controlling for number of prior mTBIs (both model terms shown):\n")
print(results, row.names = FALSE)
cat("\nAll VIFs should be well under 5 for a stable model (checked above).\n")

# Forest plot: every term in every model (predictor of interest and the
# adjustment covariate), faceted by outcome (rows) x predictor tested
# (columns). Color encodes each term's role (tested vs. controlled-for)
# rather than significance; exact p-values are printed next to each point
# instead. supp_figure_1h.R plots the same 4 underlying models with the
# tested/adjusted roles swapped, distinguished by the title.
results$term_name <- factor(results$term_name,
                             levels = c("mTBI symptom severity", "mTBI symptom frequency", "Number of prior mTBIs"))
results$p_label <- ifelse(results$p_value < 0.001, "p < 0.001", paste0("p = ", formatC(results$p_value, format = "f", digits = 3)))
results$p_label <- ifelse(results$p_value < 0.05, paste0(results$p_label, "*"), results$p_label)
# Stacked into a single column (4 rows) rather than a 2x2 grid to match the
# vertical space available in the composite manuscript figure.
results$facet_col <- paste0("mTBI ", results$predictor_label)

role_predictor <- "#2E6E9E"
role_covariate <- "#B7B7B7"
x_range <- range(c(results$ci_lower, results$ci_upper))
label_nudge <- diff(x_range) * 0.03

Forest_Plot <- ggplot(results, aes(x = estimate, y = term_name)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.12, colour = "black", linewidth = 0.6) +
  geom_point(aes(fill = term_role), shape = 21, size = 5, stroke = 0.6, colour = "black") +
  geom_text(aes(label = p_label, x = ci_upper + label_nudge), hjust = 0, size = 3.6, family = "Arial", fontface = "bold", colour = "black") +
  scale_fill_manual(values = c("Predictor of interest" = role_predictor, "Adjustment covariate" = role_covariate), name = NULL) +
  scale_y_discrete(drop = TRUE) +
  facet_grid2(outcome_label + facet_col ~ ., scales = "free_y", independent = "y") +
  coord_cartesian(xlim = c(x_range[1], x_range[2] + diff(x_range) * 0.22), clip = "off") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
    axis.line.x = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.x = element_line(colour = "black"),
    strip.text = element_text(size = 12, family = "Arial", face = "bold", colour = "black"),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
    axis.title.x = element_text(size = 12, family = "Arial", face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 11, family = "Arial", face = "bold"),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    plot.title = element_text(size = 14, family = "Arial", face = "bold", colour = "black"),
    plot.subtitle = element_text(size = 11, family = "Arial", colour = "grey30", margin = margin(b = 10)),
    plot.caption = element_text(size = 10, family = "Arial", face = "italic", colour = "black", hjust = 0, margin = margin(t = 12))
  ) +
  labs(x = "Adjusted rank-regression coefficient (95% CI)",
       title = paste(strwrap("Does mTBI symptom severity/frequency predict COVID-19 outcomes?", width = 32), collapse = "\n"),
       subtitle = paste(strwrap("Adjusted for the number of prior mTBIs", width = 40), collapse = "\n"),
       caption = "Dashed line indicates no association (coefficient = 0).")

ggsave(file.path(tif_dir, "supp_figure_1i_mtbi_severity_regression_forest.pdf"), Forest_Plot, width = 6.5, height = 13, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1i_mtbi_severity_regression_forest.pdf")), Forest_Plot, width = 6.5, height = 13, dpi = 600, device = "tiff")
