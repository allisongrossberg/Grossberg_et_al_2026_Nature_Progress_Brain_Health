# =============================================================================
# Supplementary Figure 2f – PEG vs PBS Control GFAP MFI bar graph
# =============================================================================
#
# Description: Bar plot of normalized GFAP MFI, PEG-only exosome-isolation
#   control vs. PBS negative control (rules out residual PEG from the
#   exosome isolation procedure as a confound). Welch two-sample t-test.
#
# Prerequisites: supp_figure_2fg_preprocessing_HA_TLR_Inh_MFI.R
#
# Output: supp_figure_2f_PEG_GFAP.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
source("supp_figure_2fg_preprocessing_HA_TLR_Inh_MFI.R")

library(ggplot2)
library(ggsignif)

new_labels <- c("PBS Control", "COVID-19 (-) mTBI (-) ADEs", "COVID-19 (+) mTBI (-) ADEs",
                "COVID-19 (-) mTBI (+) ADEs", "COVID-19 (+) mTBI (+) ADEs")

GFAP_result_df_clean_PEG <- GFAP_result_df_clean %>%
  filter(Condition == "PEG" | Condition == "Neg_Cont")

peg_gfap_test <- t.test(Normed_GFAP_MFI ~ Condition, data = GFAP_result_df_clean_PEG)
peg_gfap_p <- sprintf("%.4f", peg_gfap_test$p.value)

GFAP_summary_df_PEG <- GFAP_summary_df_1 %>%
  filter(Condition == "PEG" | Condition == "Neg_Cont")
GFAP_summary_df_PEG$Condition <- factor(GFAP_summary_df_PEG$Condition,
                                     levels = unique(GFAP_summary_df_PEG$Condition),
                                     labels = c("PBS Control", "PEG"))

PEG_GFAP_Plot <- ggplot(GFAP_summary_df_PEG, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "summary", width = 0.7, alpha = 0.95) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.7), width = 0.25, color = "black") +
  theme_minimal() +
  scale_fill_brewer(palette = "Greys", labels = new_labels) +
  theme(legend.title = element_blank(), legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", colour = "black"),
        axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", colour = "black"),
        legend.text = element_text(size = 18, lineheight = .9, family = "Arial", colour = "black")) +
  labs(x = "Condition", y = "Normalized GFAP MFI") +
  scale_x_discrete(labels = function(x) gsub("_", "-", x)) +
  geom_signif(
    y_position = c(max(GFAP_summary_df_PEG$Mean + GFAP_summary_df_PEG$SE) * 1.15),
    xmin = c(1), xmax = c(2),
    annotation = c(peg_gfap_p),
    tip_length = 0, color = "black", textsize = 6, size = 1
  )

print(PEG_GFAP_Plot)

path_output <- "supp_figure_2f_PEG_GFAP.pdf"
ggsave(path_output, plot = PEG_GFAP_Plot, width = 350, height = 200, units = "mm", dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", path_output), plot = PEG_GFAP_Plot, width = 350, height = 200, units = "mm", dpi = 600, device = "tiff")
