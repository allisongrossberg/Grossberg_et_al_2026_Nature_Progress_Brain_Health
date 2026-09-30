# =============================================================================
# Supplementary Figure 2b – GFAP Exosome Validation bar graph
# =============================================================================
#
# Description: Bar plot of GFAP concentration by exosome type (Simoa validation).
#   Dunn test with Holm correction for pairwise comparisons; p-values on plot.
#
# Prerequisites: supp_figure_2_input_data/GFAP_Exo_Val_Data.csv
#
# Output: supp_figure_2b_GFAP_Exo_Validation.pdf; supp_figure_2_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
GFAP_Exo_Val_Data <- read.csv("supp_figure_2_input_data/GFAP_Exo_Val_Data.csv", stringsAsFactors = FALSE)

library(dplyr)
library(dunn.test)
library(ggplot2)
library(ggsignif)

GFAP_dunn <- dunn.test(GFAP_Exo_Val_Data$Fitted_Conc_Mean, GFAP_Exo_Val_Data$Sample_Name, method = "holm")
GFAP_Exo_Validation_summary_df <- GFAP_Exo_Val_Data %>%
  group_by(Sample_Name) %>%
  summarise(
    n = n(),
    mean = mean(Fitted_Conc_Mean, na.rm = TRUE),
    se = sd(Fitted_Conc_Mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
max_value <- max(GFAP_Exo_Validation_summary_df$mean + GFAP_Exo_Validation_summary_df$se, na.rm = TRUE)
format_pvalue <- function(p) {
  if (p < 0.0001) return("p < 0.0001")
  paste("p =", formatC(p, format = "f", digits = 4))
}

# Bar plot with error bars and significance annotations (Dunn order: [1]=1vs2, [2]=1vs3, [3]=2vs3).
GFAP_Plot <- ggplot(GFAP_Exo_Validation_summary_df, aes(x = Sample_Name, y = mean, fill = Sample_Name)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.8)) +
  scale_fill_brewer(palette = "Greys") +
  theme_minimal() +
  geom_signif(
    y_position = c(max_value * 1.1, max_value * 1.2, max_value * 1.3),
    xmin = c(1, 2, 1),
    xmax = c(2, 3, 3),
    annotation = c(format_pvalue(GFAP_dunn$P.adjusted[1]), format_pvalue(GFAP_dunn$P.adjusted[3]), format_pvalue(GFAP_dunn$P.adjusted[2])),
    tip_length = 0.0, textsize = 8, size = 1
  ) +
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black")
  ) +
  labs(y = "GFAP Concentration (pg/mL)", x = "Exosome Type") +
  coord_cartesian(ylim = c(0, max_value * 1.5))

print(GFAP_Plot)

# Save derived output (summary and Dunn test results) to supp_figure_2_output_data
output_dir <- "supp_figure_2_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(GFAP_Exo_Validation_summary_df, file.path(output_dir, "supp_figure_2b_GFAP_Exo_Validation_summary.csv"), row.names = FALSE)
dunn_result_df <- data.frame(
  comparison = GFAP_dunn$comparisons,
  P_adjusted = GFAP_dunn$P.adjusted
)
write.csv(dunn_result_df, file.path(output_dir, "supp_figure_2b_GFAP_Exo_Validation_dunn.csv"), row.names = FALSE)

ggsave("supp_figure_2b_GFAP_Exo_Validation.pdf", GFAP_Plot, width = 12, height = 10, units = "in", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", "supp_figure_2b_GFAP_Exo_Validation.pdf"), GFAP_Plot, width = 12, height = 10, units = "in", dpi = 600, device = "tiff")
