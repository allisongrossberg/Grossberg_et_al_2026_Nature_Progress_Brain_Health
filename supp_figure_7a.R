# =============================================================================
# Supplementary Figure 7a - Metabolomics volcano plot, COVID-19(+)mTBI(+) vs
#   COVID-19(-)mTBI(-) ADEs (styled to match figure_7b.R's color scheme,
#   point-size-by-significance-category legend, geom_text_repel labeling,
#   theme, and naming convention)
# =============================================================================
#
# Prerequisites: run supp_figure_7a_preprocessing_2_glycolysis_reanalysis.R
#   first to produce supp_figure_7a_output_data/supp_figure_7a_combined_results.csv.
# Output: supp_figure_7a_COV_TBI_vs_Control_volcano_plot.pdf
# =============================================================================

if (basename(getwd()) == "supp_figure_7") setwd("..")
output_dir <- "supp_figure_7/supp_figure_7_output_data"
tif_dir    <- "supp_figure_7"

library(ggplot2)
library(ggrepel)
library(dplyr)

results <- read.csv(file.path(output_dir, "supp_figure_7a_combined_results.csv"), stringsAsFactors = FALSE)
target_metabolites <- c("2/3-Phospho-D-glycerate", "Phosphoenolpyruvate")

# Matches the manuscript's stated criteria for this analysis: FC > 1.0 (any
# direction) & raw p < 0.10 -- not figure_7's proteomics thresholds
# (FDR < 0.05 & |log2FC| > 2), which don't apply here (this dataset's
# fold changes top out around log2FC ~1.1; there is no magnitude gate).
results <- results %>%
  mutate(
    direction = case_when(
      p_value < 0.10 & fold_change > 1.0 ~ "Up-regulated",
      p_value < 0.10 & fold_change < 1.0 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    ),
    p_category = case_when(
      p_value < 0.001 ~ "0.001",
      p_value < 0.01  ~ "0.01",
      p_value < 0.10  ~ "0.10",
      TRUE            ~ ">=0.10"
    )
  )

# Label every significant metabolite (9 total: 7 up, 2 down) -- this
# dataset is small enough that figure_7's top-10/top-10 selection isn't
# needed; all significant hits are labeled directly.
labeled_points <- results %>% filter(direction != "Unchanged")

volcano_plot <- ggplot(results, aes(x = log2_fc, y = -log10(p_value), size = p_category)) +
  geom_point(aes(color = direction), alpha = 0.8, shape = 16) +
  scale_color_manual(values = c("Up-regulated" = "#DF4416", "Down-regulated" = "#A2B9FC", "Unchanged" = "#F1E7DB"),
                      guide = guide_legend(override.aes = list(size = 6))) +
  scale_size_manual(values = c("0.001" = 6, "0.01" = 5, "0.10" = 4, ">=0.10" = 3),
                     breaks = c("0.001", "0.01", "0.10", ">=0.10"),
                     labels = c("0.001", "0.01", "0.10", ">=0.10")) +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", color = "grey") +
  labs(title = "COV_TBI_vs_Control", x = "Log2 Fold Change", y = "-log10(p-value)",
       color = "Expression (log2FC)", size = "Significance (raw p)") +
  geom_text_repel(data = labeled_points, aes(label = metabolite), size = 5, box.padding = 1.5,
                   point.padding = 0.5, segment.color = "grey", segment.size = 0.5, max.overlaps = 25) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.title.y = element_text(size = 16, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 0, r = 30, b = 0, l = 0)),
    axis.title.x = element_text(size = 16, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 30, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(size = 16, lineheight = .9, family = "sans", colour = "black"),
    axis.text.y = element_text(size = 16, lineheight = .9, family = "sans", colour = "black"),
    legend.title = element_text(size = 16, family = "sans", face = "bold"),
    legend.text = element_text(size = 14, family = "sans"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", linewidth = 0.5, linetype = "solid")
  )

ggsave(file.path(tif_dir, "supp_figure_7a_COV_TBI_vs_Control_volcano_plot.pdf"), plot = volcano_plot,
       width = 20, height = 12, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_7a_COV_TBI_vs_Control_volcano_plot.pdf")), plot = volcano_plot,
       width = 20, height = 12, dpi = 600, device = "tiff")
cat("Saved supp_figure_7a_COV_TBI_vs_Control_volcano_plot.pdf\n")
