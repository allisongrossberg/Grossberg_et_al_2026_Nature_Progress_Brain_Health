# =============================================================================
# Supplementary Figure 1b – Total symptom frequency: combined groups
# =============================================================================
#
# Description: Box plot comparing total neurological/psychological symptom
#   frequency between "COVID-19 (+) mTBI (+)" and all other groups combined.
#   Wilcoxon rank-sum test; p-value shown on plot.
#
# Prerequisites: supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: supp_figure_1b_Neuro_Psych_Freq_Comb.pdf, supp_figure_1_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
library(dplyr)
csv_path <- "supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
if (!file.exists(csv_path)) stop("Need COAST_Study_Data_Clean_Age_Groups_add_dates.csv in supp_figure_1_input_data/.")
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(csv_path, stringsAsFactors = FALSE)

# Subset columns and rows used for severity/frequency analyses (four study groups).
Sev_Freq_Figures <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::select(participant_id, recruitment_site, age_years, qq_biological_sex, qq_group,
                recode_sev_score_total_covid, recode_sev_score_total_tbi, recode_overall_sev_score,
                covid_symptom_sev_average, tbi_symptom_sev_average, total_symptom_sev_average,
                recode_freq_score_total_covid, recode_freq_score_total_tbi, recode_overall_freq_score,
                covid_symptom_freq_average, tbi_symptom_freq_average, total_symptom_freq_average,
                qq_phq8_total_score, qq_eq5d_index_score, qq_fss_average_score, qq_fss_total_score,
                qq_gad7_total_score, qq_tbi_num,
                chronic_covid, chronic_tbi, chronic_acute_overall,
                recode_tbi_neuro_psych_sev_score, recode_tbi_neuro_psych_freq_score,
                recode_covid_neuro_psych_sev_score, recode_covid_neuro_psych_freq_score,
                recode_overall_neuro_psych_sev_score, recode_overall_neuro_psych_freq_score,
                covid_neuro_psych_symptom_sev_average, tbi_neuro_psych_symptom_sev_average,
                total_neuro_psych_symptom_sev_average, total_neuro_psych_symptom_freq_average,
                tbi_neuro_psych_symptom_freq_average, covid_neuro_psych_symptom_freq_average) %>%
  dplyr::filter(qq_group %in% c("COVID-19 (-) mTBI (-)", "COVID-19 (-) mTBI (+)",
                                "COVID-19 (+) mTBI (-)", "COVID-19 (+) mTBI (+)"))

library(ggplot2)
library(ggsignif)
library(dplyr)
library(plotrix)

# Format p-values for plot annotation (4 significant figures).
format_pvalue <- function(p) {
  if (p < 0.0001) {
    return("p < 0.0001")
  } else {
    return(paste("p =", formatC(p, format = "f", digits = 4)))
  }
}

# -----------------------------------------------------------------------------
# Prepare data: COVID-19 (+) mTBI (+) vs other groups combined
# -----------------------------------------------------------------------------

combined_data <- Sev_Freq_Figures %>%
  mutate(combined_group = factor(ifelse(qq_group == "COVID-19 (+) mTBI (+)", 
                                        "COVID-19 (+) mTBI (+)", 
                                        "Other Groups"),
                                 levels = c("Other Groups", "COVID-19 (+) mTBI (+)")))

max_value <- max(combined_data$recode_overall_neuro_psych_freq_score, na.rm = TRUE)

# Build box plot with jittered points and significance bracket.
Neuro_Psych_Freq_Comb_Final_Figure <- ggplot(combined_data, aes(x = combined_group, y = recode_overall_neuro_psych_freq_score, fill = combined_group)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(fill = combined_group), width = 0.2, shape = 21, size = 4, stroke = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Other Groups" = "#CCCCCC", "COVID-19 (+) mTBI (+)" = "#666666")) +
  theme_minimal() +

  geom_signif(y_position = max_value * 1.1,
              xmin = 1, xmax = 2,
              annotation = format_pvalue(8.131e-07),
              tip_length = 0.0, textsize = 8, size = 1) +

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
    plot.caption = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", margin = margin(t = 20), colour = "black")
  ) +
  
  labs(
    y = "Total Neurological/Psychological\nSymptom Frequency Score",
    x = "Group",
    caption = "W = 469, p < 0.0010"
  ) +
  scale_x_discrete(labels = c("Other Groups", "COVID-19 (+)\nmTBI (+)")) +
  # Extended to cover the full visible range (0 to max_value*1.2) so labels
  # reach up near this panel's single significance bracket. Shares the same
  # low-range values as figure_1d.R/supp_figure_1d.R (larger multiplier there
  # for their 6 stacked brackets).
  scale_y_continuous(breaks = c(0, 100, 200, 300)) +
  coord_cartesian(ylim = c(0, max_value * 1.2))

Neuro_Psych_Freq_Comb_Final_Figure

# Wilcoxon rank-sum test (reported in figure caption).
stat_test <- wilcox.test(recode_overall_neuro_psych_freq_score ~ combined_group, data = combined_data)
stat_test

aggregate(recode_overall_neuro_psych_freq_score ~ combined_group, 
          data = combined_data, 
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = std.error(x)))

# Export figure and plot data.
dir.create("supp_figure_1_output_data", showWarnings = FALSE)
write.csv(combined_data %>% dplyr::select(participant_id, combined_group, recode_overall_neuro_psych_freq_score), "supp_figure_1_output_data/supp_figure_1b_Neuro_Psych_Freq_Comb_data.csv", row.names = FALSE)
ggsave("supp_figure_1b_Neuro_Psych_Freq_Comb.pdf", Neuro_Psych_Freq_Comb_Final_Figure, width = 12, height = 10, units = "in", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", "supp_figure_1b_Neuro_Psych_Freq_Comb.pdf"), Neuro_Psych_Freq_Comb_Final_Figure, width = 12, height = 10, units = "in", dpi = 600, device = "tiff")
