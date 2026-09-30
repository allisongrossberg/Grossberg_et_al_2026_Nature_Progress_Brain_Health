# =============================================================================
# Figure 1g – Health-related quality of life (HRQOL, EQ-5D-5L) by group
# =============================================================================
#
# Description: Box plot of EQ-5D-5L index score across four groups.
#   Kruskal–Wallis with Dunn post-hoc (Holm); p-values on plot.
#
# Prerequisites: figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: figure_1g_EQ5D.pdf
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
) %>%
  dplyr::filter(qq_group %in% c("COVID-19 (-) mTBI (-)", "COVID-19 (-) mTBI (+)",
                                "COVID-19 (+) mTBI (-)", "COVID-19 (+) mTBI (+)"))

EQ5D_plot_data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% dplyr::select(participant_id, qq_group, qq_eq5d_index_score)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(EQ5D_plot_data %>% dplyr::select(participant_id, qq_group, qq_eq5d_index_score), file.path(output_dir, "figure_1g_EQ5D_plot_data.csv"), row.names = FALSE)

########

#EQ-5D-5L 
max_value <- max(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_eq5d_index_score, na.rm = TRUE)

# Function to format p-values to 4 significant figures
format_pvalue <- function(p) {
  if (p < 0.0001) {
    return("p < 0.0001")
  } else {
    return(paste("p =", formatC(p, format = "f", digits = 4)))
  }
}


kruskal.test(qq_eq5d_index_score ~ qq_group, data = COAST_Study_Data_Clean_Age_Groups_add_dates)
eq_dunn <- dunn.test(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_eq5d_index_score, COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group, method="holm")

aggregate(qq_eq5d_index_score ~ qq_group, 
          data = COAST_Study_Data_Clean_Age_Groups_add_dates, 
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = std.error(x)))

COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group <- factor(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group, levels = c("COVID-19 (-) mTBI (-)", "COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))
EQ_5D_5L_Figure <- ggplot(COAST_Study_Data_Clean_Age_Groups_add_dates, aes(x=qq_group, y=qq_eq5d_index_score))+
  geom_boxplot(width = 0.5, position = position_dodge(0.8), fill = "white", color = "black", alpha = 0.7, outlier.shape = NA, linewidth = 0.8) + 
  geom_jitter(aes(fill = qq_group), width = 0.2, height = 0, shape = 21, size = 4, stroke = 0.5) +
  theme_minimal()+
  scale_fill_brewer(palette = "Greys") +
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
  labs(y = "HRQOL (EQ-5D-5L)", x = "Group")+ 
  scale_x_discrete(labels = c("COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)", 
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)")) +
  coord_cartesian(ylim = c(0, max_value * 2.1)) +
  geom_signif(
    y_position = max_value * 1.15,
    xmin = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (-) mTBI (-)"),
    xmax = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (+) mTBI (-)"),
    annotation = format_pvalue(eq_dunn$P.adjusted[2]),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = max_value * 1.30,
    xmin = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (+) mTBI (-)"),
    xmax = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (-) mTBI (+)"),
    annotation = format_pvalue(eq_dunn$P.adjusted[3]),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = max_value * 1.45,
    xmin = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (-) mTBI (+)"),
    xmax = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (+) mTBI (+)"),
    annotation = format_pvalue(eq_dunn$P.adjusted[5]),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = max_value * 1.60,
    xmin = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (-) mTBI (-)"),
    xmax = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (-) mTBI (+)"),
    annotation = format_pvalue(eq_dunn$P.adjusted[1]),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = max_value * 1.75,
    xmin = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (+) mTBI (-)"),
    xmax = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (+) mTBI (+)"),
    annotation = format_pvalue(eq_dunn$P.adjusted[6]),
    tip_length = 0, color = "black", textsize = 8, size = 1
  ) +
  geom_signif(
    y_position = max_value * 1.90,
    xmin = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (-) mTBI (-)"),
    xmax = which(levels(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group) == "COVID-19 (+) mTBI (+)"),
    annotation = format_pvalue(eq_dunn$P.adjusted[4]),
    tip_length = 0, color = "black", textsize = 8, size = 1
  )



ggsave(file.path(tif_dir, "figure_1g_EQ5D.pdf"), EQ_5D_5L_Figure, width = 12, height = 10, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_1g_EQ5D.pdf")), EQ_5D_5L_Figure, width = 12, height = 10, dpi = 600, device = "tiff")
