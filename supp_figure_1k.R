# =============================================================================
# Supplementary Figure 1k – Total symptom frequency: 24 participants, 4 groups
# =============================================================================
#
# Description: Box plot of total neurological/psychological symptom frequency
#   across four study groups for a subset of 24 participants. Kruskal–Wallis
#   with Dunn post-hoc (Holm); pairwise p-values on plot.
#
# Prerequisites: supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: supp_figure_1k_Neuro_Psych_Freq_24P.pdf, supp_figure_1_output_data/*.csv
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
library(dunn.test)
library(dplyr)
library(plotrix)

# Function to format p-values to 4 significant figures
format_pvalue <- function(p) {
  if (p < 0.0001) {
    return("p < 0.0001")
  } else {
    return(paste("p =", formatC(p, format = "f", digits = 4)))
  }
}

# -----------------------------------------------------------------------------
# Subset to 24 participants and run tests
# -----------------------------------------------------------------------------

# Participant IDs for this figure (balanced representation across groups).
selected_ids <- c("84HZC", "8B2PC", "IPH6E", "PRO10", "QWMG2", "VXY1B", "XQANO", "95KXN", "B9V50", "DTMBS", "ELFHW", "P89EN", "PKXZL", "UWK9R", "9ER2B", "IO6KM", "KJHE7", "WD32N", "BW3AP", "KBPYT", "BL52P", "74U0H", "U7DJQ", "A0BND")

# Filter the dataset to include only the selected participants
selected_data <- Sev_Freq_Figures %>%
  filter(participant_id %in% selected_ids)

max_value <- max(selected_data$recode_overall_neuro_psych_freq_score, na.rm = TRUE)

kruskal.test(recode_overall_neuro_psych_freq_score ~ qq_group, data = selected_data)
neuro_psych_freq_dunn <- dunn.test(selected_data$recode_overall_neuro_psych_freq_score, selected_data$qq_group, method="holm")

aggregate(recode_overall_neuro_psych_freq_score ~ qq_group, 
          data = selected_data, 
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = std.error(x)))

selected_data$qq_group <- factor(selected_data$qq_group, levels = c("COVID-19 (-) mTBI (-)", "COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))
Neuro_Psych_Freq_24P_Final_Figure <- ggplot(selected_data, aes(x = qq_group, y = recode_overall_neuro_psych_freq_score, fill = qq_group)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(fill = qq_group), width = 0.2, shape = 21, size = 4, stroke = 0.5, alpha = 0.7) +
  scale_fill_brewer(palette = "Greys") +
  theme_minimal() +
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
    legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black")
  ) +
  
  labs(
    y = "Total Neurological/Psychological\nSymptom Frequency Score",
    x = "Group",
  ) +
  scale_x_discrete(labels = c("COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)", 
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)")) +
  # Extended to cover the full visible range (0 to max_value*2.1) so labels
  # reach up near the stacked significance brackets. Identical to figure_1d.R
  # (same max_value + multiplier).
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500)) +
  coord_cartesian(ylim = c(0, max_value * 2.1)) +
  geom_signif(
    y_position = c(max_value * 1.15),
    xmin = c(which(levels(selected_data$qq_group) == "COVID-19 (-) mTBI (-)")),
    xmax = c(which(levels(selected_data$qq_group) == "COVID-19 (+) mTBI (-)")),
    annotation = c(format_pvalue(neuro_psych_freq_dunn$P.adjusted[2])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.30),
    xmin = c(which(levels(selected_data$qq_group) == "COVID-19 (+) mTBI (-)")),
    xmax = c(which(levels(selected_data$qq_group) == "COVID-19 (-) mTBI (+)")),
    annotation = c(format_pvalue(neuro_psych_freq_dunn$P.adjusted[3])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.45),
    xmin = c(which(levels(selected_data$qq_group) == "COVID-19 (-) mTBI (+)")),
    xmax = c(which(levels(selected_data$qq_group) == "COVID-19 (+) mTBI (+)")),
    annotation = c(format_pvalue(neuro_psych_freq_dunn$P.adjusted[5])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.60),
    xmin = c(which(levels(selected_data$qq_group) == "COVID-19 (-) mTBI (-)")),
    xmax = c(which(levels(selected_data$qq_group) == "COVID-19 (-) mTBI (+)")),
    annotation = c(format_pvalue(neuro_psych_freq_dunn$P.adjusted[1])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.75),
    xmin = c(which(levels(selected_data$qq_group) == "COVID-19 (+) mTBI (-)")),
    xmax = c(which(levels(selected_data$qq_group) == "COVID-19 (+) mTBI (+)")),
    annotation = c(format_pvalue(neuro_psych_freq_dunn$P.adjusted[6])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.90),
    xmin = c(which(levels(selected_data$qq_group) == "COVID-19 (-) mTBI (-)")),
    xmax = c(which(levels(selected_data$qq_group) == "COVID-19 (+) mTBI (+)")),
    annotation = c(format_pvalue(neuro_psych_freq_dunn$P.adjusted[4])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  )

print(Neuro_Psych_Freq_24P_Final_Figure)

dir.create("supp_figure_1_output_data", showWarnings = FALSE)
write.csv(selected_data %>% dplyr::select(participant_id, qq_group, recode_overall_neuro_psych_freq_score), "supp_figure_1_output_data/supp_figure_1k_Neuro_Psych_Freq_24P_data.csv", row.names = FALSE)
ggsave("supp_figure_1k_Neuro_Psych_Freq_24P.pdf", Neuro_Psych_Freq_24P_Final_Figure, width = 12, height = 10, units = "in", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", "supp_figure_1k_Neuro_Psych_Freq_24P.pdf"), Neuro_Psych_Freq_24P_Final_Figure, width = 12, height = 10, units = "in", dpi = 600, device = "tiff")
