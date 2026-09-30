# =============================================================================
# Supplementary Figure 7d - miRNA volcano plots (all 6 ADE-group contrasts;
#   published panel uses COV_TBI_vs_Control)
#
# Description: Generates volcano plots (log2FC vs -log10 FDR) for each
#   edgeR contrast, styled identically to the proteomics volcano plots
#   (figure_7b.R): |logFC| > 2 and FDR < 0.05 define up/down-regulated
#   points, top 10 up and top 10 down labeled by significance x magnitude.
#   One TIFF per contrast.
#
#   Column names differ from the proteomics DEA output (edgeR's topTags uses
#   "FDR", not "adj.P.Val", and "miRNA" instead of "Entry_Name"), so the plot
#   function is adapted accordingly, but the visual style, thresholds, and
#   color scheme are unchanged from the proteomics template.
#
# Prerequisites: supp_figure_7_output_data/supp_figure_7d_edgeR_results_*.csv (preprocessing step 2)
#
# Output: supp_figure_7d_<contrast>_volcano_plot.pdf (one per contrast)
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "supp_figure_7") {
  if (file.exists("supp_figure_7")) setwd("supp_figure_7")
  else if (file.exists(file.path("..", "supp_figure_7"))) setwd(file.path("..", "supp_figure_7"))
}

OUTPUT_DIR <- "supp_figure_7_output_data"
edgeR_files <- list.files(OUTPUT_DIR, pattern = "^supp_figure_7d_edgeR_results_.*\\.csv$", full.names = TRUE)
if (length(edgeR_files) == 0) stop("No supp_figure_7d_edgeR_results_*.csv in ", OUTPUT_DIR, ". Run supp_figure_7d_preprocessing_2_DEA.R first.")
results_list <- lapply(edgeR_files, read.csv, stringsAsFactors = FALSE)
names(results_list) <- gsub("^supp_figure_7d_edgeR_results_|\\.csv$", "", basename(edgeR_files))

library(ggplot2)
library(ggrepel)
library(dplyr)

# Volcano plot with top up- and down-regulated miRNAs labeled (by FDR and log2FC).
# Same styling/thresholds as the proteomics volcano plot template, adapted for
# edgeR's "FDR"/"miRNA" columns instead of limma's "adj.P.Val"/"Entry_Name".
create_volcano_plot_miRNA <- function(results, title) {
  if (!("logFC" %in% colnames(results)) || !("FDR" %in% colnames(results)) || !("miRNA" %in% colnames(results))) {
    stop("The results data frame must contain 'logFC', 'FDR', and 'miRNA' columns.")
  }

  top_upregulated <- results %>%
    filter(FDR < 0.05 & logFC > 2) %>%
    mutate(score = -log10(FDR) * logFC) %>%
    slice_max(score, n = 10, with_ties = FALSE)

  top_downregulated <- results %>%
    filter(FDR < 0.05 & logFC < -2) %>%
    mutate(score = -log10(FDR) * abs(logFC)) %>%
    slice_max(score, n = 10, with_ties = FALSE)

  top_mirnas <- bind_rows(top_upregulated, top_downregulated)

  results <- results %>%
    mutate(p_value = case_when(
      FDR < 0.001 ~ "0.001",
      FDR < 0.01 & FDR >= 0.001 ~ "0.01",
      FDR < 0.05 & FDR >= 0.01 ~ "0.05",
      FDR >= 0.05 ~ ">=0.05"
    ))

  ggplot(results, aes(x = logFC, y = -log10(FDR), size = p_value)) +
    geom_point(aes(color = ifelse(FDR < 0.05 & abs(logFC) > 2,
                                  ifelse(logFC > 2, "Up-regulated", "Down-regulated"),
                                  "Unchanged")),
               alpha = 0.8, shape = 16) +
    scale_color_manual(values = c("Up-regulated" = "#DF4416", "Down-regulated" = "#A2B9FC", "Unchanged" = "#F1E7DB"),
                       guide = guide_legend(override.aes = list(size = 6))) +
    scale_size_manual(values = c("0.001" = 6, "0.01" = 5, "0.05" = 4, ">=0.05" = 3),
                      breaks = c("0.001", "0.01", "0.05", ">=0.05"),
                      labels = c("0.001", "0.01", "0.05", ">=0.05")) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(FDR)", color = "Expression (log2FC)", size = "Significance (FDR)") +
    # Extra headroom above the highest point so ggrepel has vertical room to
    # place labels for hits clustered at similar FDR (see repel note below).
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    # Significant miRNAs often cluster at nearly identical -log10(FDR), leaving
    # little vertical room to separate labels -- stronger repulsion (higher
    # force, more iterations/time) and unrestricted label placement (ylim/xlim
    # = c(NA, NA), the default clips to the panel) let ggrepel actually resolve
    # that crowding instead of leaving overlapping/crossing labels.
    geom_text_repel(data = top_mirnas, aes(label = miRNA), size = 5,
                     box.padding = 2, point.padding = 0.6,
                     segment.color = 'grey40', segment.size = 0.4,
                     min.segment.length = 0, max.overlaps = Inf,
                     force = 4, force_pull = 0.5,
                     max.iter = 100000, max.time = 5,
                     seed = 42) +
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
      legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid")
    )
}

volcano_plots <- lapply(names(results_list), function(contrast) {
  create_volcano_plot_miRNA(results_list[[contrast]], contrast)
})
names(volcano_plots) <- names(results_list)

dir.create(OUTPUT_DIR, showWarnings = FALSE)
for (contrast_name in names(volcano_plots)) {
  out_path <- paste0("supp_figure_7d_", contrast_name, "_volcano_plot.pdf")
  ggsave(out_path, plot = volcano_plots[[contrast_name]], width = 20, height = 12, dpi = 600)
  cat("Saved:", out_path, "\n")
}

cat("Volcano plots completed for", length(volcano_plots), "contrasts.\n")
