# =============================================================================
# Figure 7b – COV_TBI_vs_Control volcano plot
# =============================================================================
#
# Description: Volcano plot (log2FC vs -log10 FDR) for COV_TBI vs Control
#   limma contrast. Top proteins labeled; requires preprocessing steps 1–2.
#
# Prerequisites: figure_7_output_data (run preprocessing scripts first)
#
# Output: figure_7b_COV_TBI_vs_Control_volcano_plot.pdf
# =============================================================================
cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}
OUTPUT_DIR <- "figure_7_output_data"
results_file <- file.path(OUTPUT_DIR, "pro_HA_limma_results_COV_TBI_vs_Control.csv")
if (!file.exists(results_file)) {
  stop("Run preprocessing (steps 1–2) first. Expected: ", results_file)
}
results <- read.csv(results_file, stringsAsFactors = FALSE)

library(ggplot2)
library(ggrepel)
library(dplyr)

# Volcano plot with top up- and down-regulated proteins labeled (by FDR and log2FC).
create_volcano_plot <- function(results, title) {
  if (!("logFC" %in% colnames(results)) || !("adj.P.Val" %in% colnames(results)) || !("Entry_Name" %in% colnames(results))) {
    stop("The results data frame must contain 'logFC', 'adj.P.Val', and 'Entry_Name' columns.")
  }
  
  # Filter top 10 upregulated and downregulated by combined score (significance × magnitude)
  # Score = -log10(adj.P.Val) * abs(logFC) — higher when both significant and large fold change
  top_upregulated <- results %>%
    filter(adj.P.Val < 0.05 & logFC > 2) %>%
    mutate(score = -log10(adj.P.Val) * logFC) %>%
    slice_max(score, n = 10, with_ties = FALSE)
  
  top_downregulated <- results %>%
    filter(adj.P.Val < 0.05 & logFC < -2) %>%
    mutate(score = -log10(adj.P.Val) * abs(logFC)) %>%
    slice_max(score, n = 10, with_ties = FALSE)
  
  top_proteins <- bind_rows(top_upregulated, top_downregulated)
  
  # Map adj.P.Val to p-value categories for size legend
  results <- results %>%
    mutate(p_value = case_when(
      adj.P.Val < 0.001 ~ "0.001",
      adj.P.Val < 0.01 & adj.P.Val >= 0.001 ~ "0.01",
      adj.P.Val < 0.05 & adj.P.Val >= 0.01 ~ "0.05",
      adj.P.Val >= 0.05 ~ ">=0.05"
    ))
  
  # Create volcano plot
  ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), size = p_value)) +
    geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) > 2, 
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
    geom_text_repel(data = top_proteins, aes(label = Entry_Name), size = 5, box.padding = 1.5, point.padding = 0.5, segment.color = 'grey', segment.size = 0.5, max.overlaps = 25) +
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

vp_2 <- create_volcano_plot(results, "COV_TBI_vs_Control")
print(vp_2)

dir.create(OUTPUT_DIR, showWarnings = FALSE)
write.csv(results, file.path(OUTPUT_DIR, "figure_7b_volcano_plot_data.csv"), row.names = FALSE)

path_output <- "figure_7b_COV_TBI_vs_Control_volcano_plot.pdf"
ggsave(path_output, plot = vp_2, width = 20, height = 12, dpi = 600)
