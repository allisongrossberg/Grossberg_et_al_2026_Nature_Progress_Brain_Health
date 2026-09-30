# =============================================================================
# Supplementary Figure 1d – Anti-Aβ 1-42 IgG/IgA: combined groups
# =============================================================================
#
# Description: Box plot of anti-Aβ1-42 IgG/IgA levels comparing "COVID-19 (+)
#   mTBI (+)" to other groups combined.
#
#   Statistics: one-sided Boschloo's exact test for the pre-specified
#   hypothesis that the COVID-19 (+) mTBI (+) group has a higher proportion
#   above the control-group median than the other groups combined. Effect
#   size (odds ratio, 95% CI) comes from a two-sided Fisher's exact test on
#   the same table, since boschloo()'s own conf.int=TRUE is computationally
#   impractical at this sample size.
#
# Prerequisites: supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: supp_figure_1d_Anti_Abeta_IgG_IgA_Comb.pdf, supp_figure_1_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
csv_path <- "supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
if (!file.exists(csv_path)) stop("Need COAST_Study_Data_Clean_Age_Groups_add_dates.csv in supp_figure_1_input_data/.")
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(csv_path, stringsAsFactors = FALSE)

library(ggplot2)
library(ggsignif)
library(dplyr)
library(plotrix)
library(exact2x2)

# Subset to IgG/IgA variable and four study groups; drop missing values.
Vib_Figures_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::select(qq_group, participant_id, vib_anti_abeta_1_42_igg_iga) %>%
  filter(qq_group %in% c("COVID-19 (-) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (-)", "COVID-19 (+) mTBI (+)")) %>%
  filter(!is.na(vib_anti_abeta_1_42_igg_iga))

# Define combined-group factor and classify above/below control median.
combined_data <- Vib_Figures_Data %>%
  mutate(combined_group = factor(ifelse(qq_group == "COVID-19 (+) mTBI (+)",
                                        "COVID-19 (+) mTBI (+)",
                                        "Other Groups"),
                                 levels = c("Other Groups", "COVID-19 (+) mTBI (+)")))

# Median of control group (COVID-19 (-) mTBI (-)) for dichotomizing response.
median_value_control <- Vib_Figures_Data %>%
  filter(qq_group == "COVID-19 (-) mTBI (-)") %>%
  summarize(median_value = median(vib_anti_abeta_1_42_igg_iga, na.rm = TRUE)) %>%
  pull(median_value)
combined_data$above_median <- combined_data$vib_anti_abeta_1_42_igg_iga > median_value_control

# One-sided test: does the COVID-19 (+) mTBI (+) group have a higher
# proportion above the control median than the other groups combined?
# (pre-specified directional hypothesis; see header note.)
x_focal <- sum(combined_data$combined_group == "COVID-19 (+) mTBI (+)" & combined_data$above_median)
n_focal <- sum(combined_data$combined_group == "COVID-19 (+) mTBI (+)")
x_other <- sum(combined_data$combined_group == "Other Groups" & combined_data$above_median)
n_other <- sum(combined_data$combined_group == "Other Groups")
rate_focal <- x_focal / n_focal; rate_other <- x_other / n_other
cat(sprintf("Double+ above median: %d/%d (%.1f%%)   Other above median: %d/%d (%.1f%%)\n",
            x_focal, n_focal, 100*rate_focal, x_other, n_other, 100*rate_other))
stopifnot(rate_focal > rate_other)  # confirms direction before interpreting "greater"/"less" below

boschloo_1sided <- boschloo(x_focal, n_focal, x_other, n_other, alternative = "less")
# boschloo()'s "alternative" refers to the ratio p_other(1-p_focal)/[p_focal(1-p_other)];
# "less" (ratio < 1) corresponds to p_focal > p_other, i.e. Double+ higher - verified
# empirically above (rate_focal > rate_other) rather than assumed from the argument name.

# Effect size (OR + 95% CI) to accompany the Boschloo p-value; boschloo()'s
# own conf.int=TRUE is impractical at this n.
fisher_2sided <- fisher.test(
  matrix(c(x_focal, n_focal - x_focal, x_other, n_other - x_other), nrow = 2, byrow = TRUE,
         dimnames = list(c("Double+", "Other"), c("above", "below")))
)
or_estimate <- unname(fisher_2sided$estimate)
or_ci_lower <- fisher_2sided$conf.int[1]
or_ci_upper <- fisher_2sided$conf.int[2]

cat(sprintf("One-sided Boschloo's exact (Double+ > Other): p = %.4f\n", boschloo_1sided$p.value))
cat(sprintf("Two-sided Fisher's OR (effect size for caption): %.2f, 95%% CI [%.2f, %.2f]\n",
            or_estimate, or_ci_lower, or_ci_upper))

format_pvalue <- function(p) {
  if (p < 0.0001) return("p < 0.0001")
  return(paste("p =", formatC(p, format = "f", digits = 4)))
}
format_ci <- function(lower, upper) {
  upper_str <- if (is.infinite(upper)) "∞" else sprintf("%.2f", upper)
  sprintf("[%.2f, %s]", lower, upper_str)
}

max_value <- max(combined_data$vib_anti_abeta_1_42_igg_iga, na.rm = TRUE)

# Wrapped to two lines: the full caption is wider than the 12in canvas at
# this font size regardless of alignment.
caption_text <- sprintf("Boschloo's exact (one-sided): %d/%d vs %d/%d above median;\nOR = %.2f, 95%% CI %s; %s",
                         x_focal, n_focal, x_other, n_other, or_estimate, format_ci(or_ci_lower, or_ci_upper),
                         format_pvalue(boschloo_1sided$p.value))

Anti_Abeta_IgG_IgA_Comb_Figure <- ggplot(combined_data, aes(x = combined_group, y = vib_anti_abeta_1_42_igg_iga, fill = combined_group)) +
  # Control-group median: the threshold the Boschloo's exact test above is computed against.
  geom_hline(yintercept = median_value_control, linetype = "dotted", color = "#DD4726", linewidth = 3) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(fill = combined_group), width = 0.2, shape = 21, size = 4, stroke = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Other Groups" = "#CCCCCC", "COVID-19 (+) mTBI (+)" = "#666666")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    # Left-aligned: this caption is too wide for a right-anchored (hjust=1)
    # caption without overflowing the left edge.
    plot.caption = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", margin = margin(t = 20), colour = "black", hjust = 0)
  ) +
  labs(
    y = "Anti-Aβ1-42 (IgG/IgA)",
    x = "Group",
    caption = caption_text
  ) +
  scale_x_discrete(labels = c("Other Groups", "COVID-19 (+)\nmTBI (+)")) +
  # Explicit breaks so this axis reads identically to figure_1j.R (same
  # variable, same pooled data).
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
  coord_cartesian(ylim = c(0, max_value * 1.2))

print(Anti_Abeta_IgG_IgA_Comb_Figure)

dir.create("supp_figure_1_output_data", showWarnings = FALSE)
write.csv(combined_data %>% dplyr::select(participant_id, combined_group, vib_anti_abeta_1_42_igg_iga), "supp_figure_1_output_data/supp_figure_1d_Anti_Abeta_IgG_IgA_Comb_data.csv", row.names = FALSE)
stats_out <- data.frame(
  focal_group = "COVID-19 (+) mTBI (+)", comparison_group = "Other 3 groups pooled",
  threshold_type = "control-group median", threshold_value = median_value_control,
  x_focal = x_focal, n_focal = n_focal, x_other = x_other, n_other = n_other,
  boschloo_p_onesided = boschloo_1sided$p.value,
  or_twosided = or_estimate, or_ci_lower_twosided = or_ci_lower, or_ci_upper_twosided = or_ci_upper
)
write.csv(stats_out, "supp_figure_1_output_data/supp_figure_1d_stats.csv", row.names = FALSE)
ggsave("supp_figure_1d_Anti_Abeta_IgG_IgA_Comb.pdf", Anti_Abeta_IgG_IgA_Comb_Figure, width = 12, height = 10, units = "in", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", "supp_figure_1d_Anti_Abeta_IgG_IgA_Comb.pdf"), Anti_Abeta_IgG_IgA_Comb_Figure, width = 12, height = 10, units = "in", dpi = 600, device = "tiff")
