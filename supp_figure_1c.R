# =============================================================================
# Supplementary Figure 1c – Plasma Aβ42/Aβ40 ratio: combined groups
# =============================================================================
#
# Description: Box plot of plasma Aβ42/Aβ40 ratio (N3PA assay) comparing
#   "COVID-19 (+) mTBI (+)" to other groups combined. Reference line at 0.04
#   (literature-derived amyloid-abnormality cutoff).
#
#   Statistics: one-sided Boschloo's exact test for the pre-specified
#   hypothesis that the COVID-19 (+) mTBI (+) group has a higher proportion
#   below 0.04 than the other groups combined. Effect size (odds ratio,
#   95% CI) comes from a two-sided Fisher's exact test on the same table,
#   since boschloo()'s own conf.int=TRUE is computationally impractical at
#   this sample size.
#
# Prerequisites: supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv,
#   supp_figure_1_input_data/N3PA_SRX_Data_7_30_24.xlsx
#
# Output: supp_figure_1c_Abeta_Ratio_Comb.pdf, supp_figure_1_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
library(dplyr)
csv_path <- "supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
if (!file.exists(csv_path)) stop("Need COAST_Study_Data_Clean_Age_Groups_add_dates.csv in supp_figure_1_input_data/.")
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(csv_path, stringsAsFactors = FALSE)
join_and_select <- function(df_A, df_B, columns_to_select) {
  if (!("participant_id" %in% names(df_A)) || !("participant_id" %in% names(df_B))) stop("Both dataframes must have 'participant_id'")
  if (!all(columns_to_select %in% names(df_A))) stop("columns_to_select must be in df_A")
  df_A %>% dplyr::select(participant_id, dplyr::all_of(columns_to_select)) %>% dplyr::right_join(df_B, by = "participant_id")
}

library(ggplot2)
library(ggsignif)
library(dplyr)
library(plotrix)
library(readxl)
library(exact2x2)

# Load N3PA assay data and retain only complete samples.
n3pa_csv  <- "supp_figure_1_input_data/N3PA_SRX_Data_7_30_24.csv"
n3pa_xlsx <- "supp_figure_1_input_data/N3PA_SRX_Data_7_30_24.xlsx"
if (file.exists(n3pa_csv)) {
  N3PA_SRX_Data <- read.csv(n3pa_csv, stringsAsFactors = FALSE)
} else {
  N3PA_SRX_Data <- read_excel(n3pa_xlsx)
}
N3PA_SRX_Data_Clean <- N3PA_SRX_Data %>% filter((N3PA_Status == "Complete"))

N3PA_SRX_Data_Clean_join <- join_and_select(COAST_Study_Data_Clean_Age_Groups_add_dates,
                                            N3PA_SRX_Data_Clean,
                                            c("qq_biological_sex", "age_years", "qq_group",
                                              "recode_overall_neuro_psych_sev_score"))

N3PA_SRX_Data_Clean_join$Abeta42_Abeta40_Ratio_R <- N3PA_SRX_Data_Clean_join$N3PA_Abeta42_Concentration / N3PA_SRX_Data_Clean_join$N3PA_Abeta40_Concentration
N3PA_SRX_Data_Clean_join <- N3PA_SRX_Data_Clean_join %>% filter(!is.na(qq_group), !is.na(Abeta42_Abeta40_Ratio_R))

# Merge with COAST demographics and define combined-group factor.
combined_data <- N3PA_SRX_Data_Clean_join %>%
  mutate(combined_group = factor(ifelse(qq_group == "COVID-19 (+) mTBI (+)",
                                        "COVID-19 (+) mTBI (+)",
                                        "Other Groups"),
                                 levels = c("Other Groups", "COVID-19 (+) mTBI (+)")))

aggregate(Abeta42_Abeta40_Ratio_R ~ combined_group, data = combined_data,
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = std.error(x)))

# One-sided test: does the COVID-19 (+) mTBI (+) group have a higher
# proportion below 0.04 than the other groups combined? (pre-specified
# directional hypothesis; see header note.)
n_other <- sum(combined_data$combined_group == "Other Groups")
n_focal <- sum(combined_data$combined_group == "COVID-19 (+) mTBI (+)")
x_other <- sum(combined_data$combined_group == "Other Groups" & combined_data$Abeta42_Abeta40_Ratio_R < 0.04)
x_focal <- sum(combined_data$combined_group == "COVID-19 (+) mTBI (+)" & combined_data$Abeta42_Abeta40_Ratio_R < 0.04)
rate_focal <- x_focal / n_focal; rate_other <- x_other / n_other
cat(sprintf("Double+ below 0.04: %d/%d (%.1f%%)   Other below 0.04: %d/%d (%.1f%%)\n",
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
         dimnames = list(c("Double+", "Other"), c("below", "above")))
)
or_estimate <- unname(fisher_2sided$estimate)
or_ci_lower <- fisher_2sided$conf.int[1]
or_ci_upper <- fisher_2sided$conf.int[2]

cat(sprintf("One-sided Boschloo's exact (Double+ > Other): p = %.4f\n", boschloo_1sided$p.value))
cat(sprintf("Two-sided Fisher's OR (effect size for caption): %.2f, 95%% CI [%.2f, %.2f]\n",
            or_estimate, or_ci_lower, or_ci_upper))

n_abeta <- nrow(combined_data)

format_pvalue <- function(p) {
  if (p < 0.0001) return("p < 0.0001")
  return(paste("p =", formatC(p, format = "f", digits = 4)))
}
format_ci <- function(lower, upper) {
  upper_str <- if (is.infinite(upper)) "∞" else sprintf("%.2f", upper)
  sprintf("[%.2f, %s]", lower, upper_str)
}

max_value <- max(combined_data$Abeta42_Abeta40_Ratio_R, na.rm = TRUE)

# Upper bound set from the data (30% headroom, rounded up to the next 0.01)
# rather than hardcoded, so no participants are clipped from the plot.
y_upper <- ceiling(max_value * 1.3 / 0.01) * 0.01

# Wrapped to two lines: the full caption is wider than the 12in canvas at
# this font size regardless of alignment.
caption_text <- sprintf("Boschloo's exact (one-sided): %d/%d vs %d/%d below 0.04;\nOR = %.2f, 95%% CI %s; %s",
                         x_focal, n_focal, x_other, n_other, or_estimate, format_ci(or_ci_lower, or_ci_upper),
                         format_pvalue(boschloo_1sided$p.value))

Abeta_Ratio_Comb_Figure <- ggplot(combined_data, aes(x = combined_group, y = Abeta42_Abeta40_Ratio_R, fill = combined_group)) +
  geom_hline(yintercept = 0.04, linetype = "dotted", color = "#DD4726", linewidth = 3) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(fill = combined_group), width = 0.2, shape = 21, size = 4, stroke = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Other Groups" = "#CCCCCC", "COVID-19 (+) mTBI (+)" = "#666666")) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, y_upper), breaks = seq(0, y_upper, by = 0.01)) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    # Left-aligned to match supp_figure_1d's caption.
    plot.caption = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", margin = margin(t = 20), colour = "black", hjust = 0)
  ) +
  labs(
    y = "Plasma Aβ42/Aβ40 Ratio",
    x = "Group",
    caption = caption_text
  ) +
  scale_x_discrete(labels = c("Other Groups", "COVID-19 (+)\nmTBI (+)"))

print(Abeta_Ratio_Comb_Figure)

dir.create("supp_figure_1_output_data", showWarnings = FALSE)
write.csv(combined_data %>% dplyr::select(participant_id, combined_group, Abeta42_Abeta40_Ratio_R), "supp_figure_1_output_data/supp_figure_1c_Abeta_Ratio_Comb_data.csv", row.names = FALSE)
stats_out <- data.frame(
  focal_group = "COVID-19 (+) mTBI (+)", comparison_group = "Other 3 groups pooled",
  threshold = 0.04, x_focal = x_focal, n_focal = n_focal, x_other = x_other, n_other = n_other,
  boschloo_p_onesided = boschloo_1sided$p.value,
  or_twosided = or_estimate, or_ci_lower_twosided = or_ci_lower, or_ci_upper_twosided = or_ci_upper
)
write.csv(stats_out, "supp_figure_1_output_data/supp_figure_1c_stats.csv", row.names = FALSE)
ggsave("supp_figure_1c_Abeta_Ratio_Comb.pdf", Abeta_Ratio_Comb_Figure, width = 12, height = 10, units = "in", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", "supp_figure_1c_Abeta_Ratio_Comb.pdf"), Abeta_Ratio_Comb_Figure, width = 12, height = 10, units = "in", dpi = 600, device = "tiff")
