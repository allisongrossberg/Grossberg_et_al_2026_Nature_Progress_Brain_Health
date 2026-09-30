# =============================================================================
# Figure 1d – Total neurological/psychological symptom frequency score by group
# =============================================================================
#
# Description: Box plot of total neuro/psych symptom frequency across four
#   groups. Kruskal–Wallis with Dunn post-hoc (Holm); p-values on plot.
#
# Prerequisites: figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: figure_1d_Neuro_Psych_Freq.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
output_dir <- if (dir.exists("figure_1")) "figure_1/figure_1_output_data" else "figure_1_output_data"
tif_dir    <- if (dir.exists("figure_1")) "figure_1" else "."
library(dplyr)
library(ggplot2)
library(ggsignif)
library(dunn.test)
library(plotrix)

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(
  file.path(output_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"),
  stringsAsFactors = FALSE
)

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

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Sev_Freq_Figures %>% dplyr::select(participant_id, qq_group, recode_overall_neuro_psych_freq_score), file.path(output_dir, "figure_1d_Sev_Freq_Figures.csv"), row.names = FALSE)

#######
#Symptom Frequency Figure 
max_value <- max(Sev_Freq_Figures$recode_overall_neuro_psych_freq_score, na.rm = TRUE)

# Function to format p-values to 4 significant figures
format_pvalue <- function(p) {
  if (p < 0.0001) {
    return("p < 0.0001")
  } else {
    return(paste("p =", formatC(p, format = "f", digits = 4)))
  }
}

kruskal.test(recode_overall_neuro_psych_freq_score ~ qq_group, data = Sev_Freq_Figures)
Freq_Figures_dunn <- dunn.test(Sev_Freq_Figures$recode_overall_neuro_psych_freq_score, Sev_Freq_Figures$qq_group, method="holm")
aggregate(Sev_Freq_Figures$recode_overall_neuro_psych_freq_score,by=list(Sev_Freq_Figures$qq_group), FUN=mean,na.rm=TRUE) #range

aggregate(recode_overall_neuro_psych_freq_score ~ qq_group, 
          data = Sev_Freq_Figures, 
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = std.error(x)))

Sev_Freq_Figures$qq_group <- factor(Sev_Freq_Figures$qq_group, levels = c("COVID-19 (-) mTBI (-)", "COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))
Neuro_Psych_Freq_Final_Figure <- ggplot(Sev_Freq_Figures, aes(x = qq_group, y = recode_overall_neuro_psych_freq_score)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.8), fill = "white", color = "black", alpha = 0.7, outlier.shape = NA, linewidth = 0.8) + 
  geom_jitter(aes(fill = qq_group), width = 0.2, height = 0, shape = 21, size = 4, stroke = 0.5) +
  scale_fill_brewer(palette = "Greys") +
  theme_minimal() +
  theme(legend.title=element_blank(),legend.position="none", 
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"), 
        axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"), 
        legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"))+ 
  labs(
    y = "Total Neurological/Psychological\nSymptom Frequency Score",
    x = "Group"
  ) +
  scale_x_discrete(labels = c("COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)",
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)")) +
  # Extended to cover the full visible range (0 to max_value*2.1) so labels
  # reach up near the stacked significance brackets. Shares the same
  # low-range values as efigure_1b.R (smaller multiplier, 1 bracket).
  # Identical to efigure_1d.R (same max_value + multiplier).
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500)) +
  coord_cartesian(ylim = c(0, max_value * 2.1)) +
  geom_signif(
    y_position = c(max_value * 1.15),
    xmin = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (-) mTBI (-)")),
    xmax = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (+) mTBI (-)")),
    annotation = c(format_pvalue(Freq_Figures_dunn$P.adjusted[2])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.30),
    xmin = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (+) mTBI (-)")),
    xmax = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (-) mTBI (+)")),
    annotation = c(format_pvalue(Freq_Figures_dunn$P.adjusted[3])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.45),
    xmin = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (-) mTBI (+)")),
    xmax = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (+) mTBI (+)")),
    annotation = c(format_pvalue(Freq_Figures_dunn$P.adjusted[5])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.60),
    xmin = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (-) mTBI (-)")),
    xmax = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (-) mTBI (+)")),
    annotation = c(format_pvalue(Freq_Figures_dunn$P.adjusted[1])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.75),
    xmin = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (+) mTBI (-)")),
    xmax = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (+) mTBI (+)")),
    annotation = c(format_pvalue(Freq_Figures_dunn$P.adjusted[6])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = c(max_value * 1.90),
    xmin = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (-) mTBI (-)")),
    xmax = c(which(levels(Sev_Freq_Figures$qq_group) == "COVID-19 (+) mTBI (+)")),
    annotation = c(format_pvalue(Freq_Figures_dunn$P.adjusted[4])),
    tip_length = 0, color = "black", textsize = 8, size = 1
  )

ggsave(file.path(tif_dir, "figure_1d_Neuro_Psych_Freq.pdf"), Neuro_Psych_Freq_Final_Figure, width = 12, height = 10, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_1d_Neuro_Psych_Freq.pdf")), Neuro_Psych_Freq_Final_Figure, width = 12, height = 10, dpi = 600, device = "tiff")
