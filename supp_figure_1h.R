# =============================================================================
# Supplementary Figure 1h – Does number of prior mTBIs predict COVID-19
#   outcomes after controlling for the severity/frequency of the injury?
# =============================================================================
#
# Description: supp_figure_1h_preprocessing_mtbi_count.R shows unadjusted
#   Wilcoxon/Spearman comparisons of Single (1) vs. Multiple (2+) prior
#   mTBIs against COVID-19 outcomes. A participant with more prior mTBIs
#   may also tend to have had a more severe/persistent one, so this adjusts
#   the count effect for the severity/frequency of the specific target mTBI
#   (tbi_severity_target_incidence / tbi_frequency_target_incidence, from
#   supp_figure_1_preprocessing_mtbi_severity_target.R). Their correlation
#   with n_tbi_before_covid is weak (rho = 0.18 severity, 0.10 frequency),
#   so this is a check on a modest confound, not an expected reversal.
#
#   Same population as supp_figure_1h_preprocessing_mtbi_count.R:
#   double-exposed, single COVID-19 infection, >=1 pre-COVID-19 mTBI
#   (n = 40). Count is used continuously (raw n_tbi_before_covid) rather
#   than re-binned, since binning both count and severity in the same
#   small model would fragment the sample further.
#
# Prerequisites: supp_figure_1/supp_figure_1_output_data/supp_figure_1_mtbi_count_input_data.csv
#   (run supp_figure_1h_preprocessing_mtbi_count.R first); supp_figure_1/supp_figure_1_input_data/mtbi_severity_target_incidence.csv
#   (run supp_figure_1_preprocessing_mtbi_severity_target.R first).
#
# Output: supp_figure_1h_mtbi_count_regression_forest.pdf; adjusted model
#   coefficients in supp_figure_1_output_data/supp_figure_1h_mtbi_count_regression_stats.csv
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

path_input <- file.path(output_dir, "supp_figure_1_mtbi_count_input_data.csv")
path_sev   <- file.path(input_dir, "mtbi_severity_target_incidence.csv")
missing <- c(
  if (!file.exists(path_input)) "supp_figure_1_mtbi_count_input_data.csv (run supp_figure_1h_preprocessing_mtbi_count.R first)",
  if (!file.exists(path_sev))   "mtbi_severity_target_incidence.csv (run supp_figure_1_preprocessing_mtbi_severity_target.R first)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(ggplot2)
library(car)  # vif()
library(ggh4x)  # facet_grid2() for per-panel y-axis dropping

Corr_Data <- read.csv(path_input, stringsAsFactors = FALSE)
Severity_Data <- read.csv(path_sev, stringsAsFactors = FALSE)

Corr_Data <- Corr_Data %>%
  dplyr::left_join(Severity_Data %>% dplyr::select(participant_id, tbi_severity_target_incidence,
                                                     tbi_frequency_target_incidence),
                    by = "participant_id")

cat("n =", nrow(Corr_Data), "; complete severity/frequency:",
    sum(stats::complete.cases(Corr_Data$tbi_severity_target_incidence, Corr_Data$tbi_frequency_target_incidence)), "\n")

cat("\nSpearman correlations among candidate predictors (checking for collinearity before adjusting):\n")
print(round(cor(Corr_Data[, c("n_tbi_before_covid", "tbi_severity_target_incidence", "tbi_frequency_target_incidence")],
                 method = "spearman", use = "pairwise.complete.obs"), 2))

outcome_specs <- list(
  list(col = "recode_covid_neuro_psych_sev_score",  label = "COVID Severity"),
  list(col = "recode_covid_neuro_psych_freq_score", label = "COVID Frequency"),
  list(col = "Abeta42_Abeta40_Ratio_R",              label = "Aβ42/Aβ40"),
  list(col = "vib_anti_abeta_1_42_igg_iga",           label = "Anti-Aβ IgG/IgA")
)
covariate_specs <- list(
  list(col = "tbi_severity_target_incidence",  label = "Severity"),
  list(col = "tbi_frequency_target_incidence", label = "Frequency")
)

# Fits rank(outcome) ~ rank(n_tbi_before_covid) + rank(covariate) and
# returns BOTH terms - count (the predictor of interest) AND the
# severity/frequency covariate - so the full model is visible.
fit_adjusted_model <- function(df, outcome_col, covariate_col, covariate_name) {
  sub <- df %>%
    dplyr::select(all_of(c("n_tbi_before_covid", outcome_col, covariate_col))) %>%
    stats::na.omit()
  sub_ranked <- as.data.frame(lapply(sub, rank))
  form <- stats::as.formula(paste0(outcome_col, " ~ n_tbi_before_covid + ", covariate_col))
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
    term_row("n_tbi_before_covid", "Number of prior mTBIs", "Predictor of interest"),
    term_row(covariate_col, covariate_name, "Adjustment covariate")
  )
}

results <- do.call(rbind, lapply(covariate_specs, function(cv) {
  do.call(rbind, lapply(outcome_specs, function(o) {
    res <- fit_adjusted_model(Corr_Data, o$col, cv$col, paste0("mTBI symptom ", tolower(cv$label)))
    res$covariate_label <- cv$label
    res$outcome_label <- o$label
    res
  }))
}))
# FDR correction applies only to the count (predictor-of-interest) tests -
# the covariate's own p-value is shown for transparency, not as a formal test.
results$p_fdr <- NA_real_
is_predictor <- results$term_role == "Predictor of interest"
results$p_fdr[is_predictor] <- p.adjust(results$p_value[is_predictor], method = "fdr")
results <- results[, c("covariate_label", "outcome_label", "term_role", "term_name", "n",
                        "estimate", "ci_lower", "ci_upper", "p_value", "p_fdr", "r_squared", "vif")]

write.csv(results, file.path(output_dir, "supp_figure_1h_mtbi_count_regression_stats.csv"), row.names = FALSE)
cat("\nAdjusted (rank-regression) effect of # of prior mTBIs on COVID/biomarker outcomes,\n",
    "controlling for severity or frequency of the target incidence (both model terms shown):\n")
print(results, row.names = FALSE)
cat("\nAll VIFs should be well under 5 for a stable model (checked above).\n")

# -----------------------------------------------------------------------
# Forest plot: EVERY term in every model (count, the predictor of interest,
# AND the severity/frequency covariate), for COVID Severity/Frequency only
# (the two outcomes with any unadjusted signal elsewhere in this folder;
# biomarkers are in the stats CSV above for completeness but not plotted).
# Faceted by outcome (rows) x covariate adjusted for (columns). Color
# encodes term ROLE (predictor of interest vs. adjustment covariate), not
# significance - exact p-values are printed next to each point instead of a
# color cutoff. Title/subtitle state which term is the tested predictor here
# (count), since this figure and supp_figure_1_mtbi_severity_regression.R's plot
# the same 4 underlying models with the tested/adjusted roles swapped - the
# only way to tell the two figures apart is to read the title.
# -----------------------------------------------------------------------
plot_data <- results %>% dplyr::filter(outcome_label %in% c("COVID Severity", "COVID Frequency"))
plot_data$term_name <- factor(plot_data$term_name,
                               levels = c("Number of prior mTBIs", "mTBI symptom severity", "mTBI symptom frequency"))
plot_data$p_label <- ifelse(plot_data$p_value < 0.001, "p < 0.001", paste0("p = ", formatC(plot_data$p_value, format = "f", digits = 3)))
plot_data$p_label <- ifelse(plot_data$p_value < 0.05, paste0(plot_data$p_label, "*"), plot_data$p_label)
# Stacked into a single column (4 rows) rather than a 2x2 grid, so the figure
# grows taller/narrower instead of wider - matches the vertical space
# available in the composite manuscript figure. Strip text is just "mTBI
# Frequency/Severity" (not "Adjusted for mTBI ...") since the "adjusted for"
# framing is already stated once in the subtitle and the longer phrasing
# doesn't fit in the narrower strip width.
plot_data$facet_col <- paste0("mTBI ", plot_data$covariate_label)

role_predictor <- "#2E6E9E"
role_covariate <- "#B7B7B7"
x_range <- range(c(plot_data$ci_lower, plot_data$ci_upper))
label_nudge <- diff(x_range) * 0.03

Forest_Plot <- ggplot(plot_data, aes(x = estimate, y = term_name)) +
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
       title = paste(strwrap("Does the number of prior mTBIs predict COVID-19 outcomes?", width = 32), collapse = "\n"),
       subtitle = paste(strwrap("Adjusted for how severe or frequent that mTBI's symptoms were", width = 40), collapse = "\n"),
       caption = "Dashed line indicates no association (coefficient = 0).")

ggsave(file.path(tif_dir, "supp_figure_1h_mtbi_count_regression_forest.pdf"), Forest_Plot, width = 6.5, height = 13, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1h_mtbi_count_regression_forest.pdf")), Forest_Plot, width = 6.5, height = 13, dpi = 600, device = "tiff")
