# =============================================================================
# Supplementary Figure 5d – mtDNA D-Loop concentration: astrocyte-derived exosomes
# =============================================================================
#
# Description: Box plot of mtDNA D-Loop concentration by study group
#   (astrocyte-derived exosomes). One-way ANOVA; caption reports F and p.
#
# Prerequisites: supp_figure_5_input_data/mtDNA_ADE_DLoop_Results_join_ADEs.csv
#
# Output: supp_figure_5d_mtDNA_D_Loop_ADE.pdf; supp_figure_5_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
mtDNA_ADE_DLoop_Results_join_ADEs <- read.csv("supp_figure_5_input_data/mtDNA_ADE_DLoop_Results_join_ADEs.csv", stringsAsFactors = FALSE)

library(ggplot2)
library(car)
library(plotrix)
library(dplyr)

# -----------------------------------------------------------------------------
# Astrocyte-derived exosome figures
# -----------------------------------------------------------------------------
mtDNA_ADE_DLoop_Results_join_ADEs$qq_group <- factor(mtDNA_ADE_DLoop_Results_join_ADEs$qq_group, levels = c("COVID-19 (-) mTBI (-)", "COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))
Figure_mtDNA_D_Loop_ADE <- ggplot(mtDNA_ADE_DLoop_Results_join_ADEs, aes(x = qq_group, y = D_Loop_Conc_copies_µL_, fill = qq_group)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, shape = 21, size = 3, stroke = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("#FFFFFF", "#D0D0D0", "#808080", "#080808")) +
  theme_minimal() +
  
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),
    panel.background = element_rect(fill = '#FFFFFF', color = '#FFFFFF'),
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
    plot.caption = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", margin = margin(t = 20), colour = "black")
  ) +
  
  labs(
    y = "D-Loop Concentration (copies/µL)",
    x = "Group",
    caption = "F(3, 21) = 0.70, p = 0.5650"
  ) +
  scale_x_discrete(labels = c("COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)", 
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)"),
                  drop = FALSE)
Figure_mtDNA_D_Loop_ADE

aov_result <- aov(D_Loop_Conc_copies_µL_ ~ qq_group, data = mtDNA_ADE_DLoop_Results_join_ADEs)
summary(aov_result)

# Shapiro-Wilk test for normality
shapiro.test(residuals(aov_result))

# Levene's test for homogeneity of variances
leveneTest(D_Loop_Conc_copies_µL_ ~ qq_group, data = mtDNA_ADE_DLoop_Results_join_ADEs)

# Summary and ANOVA output (derived data, not same as input)
output_dir <- "supp_figure_5_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
summary_df <- mtDNA_ADE_DLoop_Results_join_ADEs %>%
  group_by(qq_group) %>%
  summarise(n = n(), mean = mean(D_Loop_Conc_copies_µL_, na.rm = TRUE), se = std.error(D_Loop_Conc_copies_µL_), .groups = "drop")
s_aov <- summary(aov_result)
anova_df <- data.frame(F_value = s_aov[[1]]$`F value`[1], df1 = s_aov[[1]]$Df[1], df2 = s_aov[[1]]$Df[2], p_value = s_aov[[1]]$`Pr(>F)`[1])
write.csv(summary_df, file.path(output_dir, "supp_figure_5d_D_Loop_ADE_summary.csv"), row.names = FALSE)
write.csv(anova_df, file.path(output_dir, "supp_figure_5d_D_Loop_ADE_anova.csv"), row.names = FALSE)

ggsave("supp_figure_5d_mtDNA_D_Loop_ADE.pdf", Figure_mtDNA_D_Loop_ADE, width = 12, height = 10, units = "in", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", "supp_figure_5d_mtDNA_D_Loop_ADE.pdf"), Figure_mtDNA_D_Loop_ADE, width = 12, height = 10, units = "in", dpi = 600, device = "tiff")
