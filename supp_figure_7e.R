# =============================================================================
# Supplementary Figure 7e - miRNA-target regulatory network (COV_TBI vs Control)
# =============================================================================
#
# Description: Network linking the 7 miRNAs significant in
#   COV_TBI_vs_Control (left) to genes on the right, in two blocks: (1) the
#   54 genes multiMiR lists as a validated target of at least one of those
#   miRNAs, connected by an edge colored by directional consistency (teal =
#   miRNA/target move in opposite directions, consistent with direct
#   repression; orange = same direction); (2) the 22 top differentially
#   expressed proteins with no validated miRNA target and no
#   reactive-astrocyte marker citation, appended below with no edge, so
#   their absence of a connection is visible. Gene nodes are colored by
#   their own logFC and labeled with logFC + significance stars.
#
#   Caveat: validated-target status is evidence from other experimental
#   systems, not proof these miRNAs regulate these proteins in this
#   dataset -- report any edge as "consistent with a candidate regulatory
#   relationship," not causality.
#
# Prerequisites: supp_figure_7_output_data/supp_figure_7e_mirna_target_overlap_all.csv
#   (from supp_figure_7e_preprocessing_1_mirna_target_overlap.R);
#   supp_figure_7_output_data/supp_figure_7e_top25_significant_results.csv
#   (from supp_figure_7e_preprocessing_2_top25_significant_results.R)
#
# Output: supp_figure_7e_top25_mirna_combined.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
if (basename(getwd()) == "supp_figure_7") setwd("..")

output_dir <- "supp_figure_7/supp_figure_7_output_data"
tif_dir    <- "supp_figure_7"

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(ggnewscale)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

d <- read.csv(file.path(output_dir, "supp_figure_7e_mirna_target_overlap_all.csv"), stringsAsFactors = FALSE)

top25_path <- file.path(output_dir, "supp_figure_7e_top25_significant_results.csv")
if (!file.exists(top25_path)) stop("Run supp_figure_7e_preprocessing_2_top25_significant_results.R first. Expected: ", top25_path)
top25 <- read.csv(top25_path, stringsAsFactors = FALSE) %>%
  dplyr::transmute(gene_symbol, logFC = logFC_COV_TBI, adj.P.Val = adj.P.Val_COV_TBI)

# One edge per miRNA-target pair (drop the duplicate "Escartin panel" row for
# ICAM1 -- that pair is already present via the DE-protein-arm row).
edges <- d %>%
  group_by(miRNA, target_symbol) %>%
  summarise(
    miRNA_logFC = dplyr::first(miRNA_logFC),
    logFC_target = dplyr::first(logFC_target),
    adj.P.Val_target = dplyr::first(adj.P.Val_target),
    directionality_consistent = dplyr::first(directionality_consistent),
    .groups = "drop"
  )

# Left column: order miRNAs down-regulated first, then up-regulated
mirna_order <- edges %>% distinct(miRNA, miRNA_logFC) %>% arrange(miRNA_logFC) %>% pull(miRNA)
mirna_logfc <- edges %>% distinct(miRNA, miRNA_logFC) %>% arrange(miRNA_logFC) %>% pull(miRNA_logFC)
left_rank <- setNames(seq_along(mirna_order), mirna_order)
n_left <- length(mirna_order)

# Block 1 (validated targets): grouped by the mean left-rank of connected
# miRNA(s), then by own logFC -- keeps each miRNA's targets together and
# places shared targets sensibly in between, minimizing edge crossings.
target_key <- edges %>%
  group_by(target_symbol) %>%
  summarise(rank_key = mean(left_rank[miRNA]), logFC_target = dplyr::first(logFC_target),
            adj.P.Val_target = dplyr::first(adj.P.Val_target), .groups = "drop") %>%
  arrange(rank_key, desc(logFC_target))
target_order <- target_key$target_symbol
n_target <- length(target_order)

# Block 2 (top25-only, no validated target relationship): appended below
# block 1, ranked by combined FC x significance score (same convention as
# supp_figure_7e_preprocessing_2_top25_significant_results.R's own gene selection), descending
rank_score <- function(logfc, padj) -log10(pmax(padj, 1e-300)) * logfc
extra_order <- top25 %>%
  filter(!gene_symbol %in% target_order) %>%
  mutate(score = rank_score(logFC, adj.P.Val)) %>%
  arrange(desc(score)) %>%
  pull(gene_symbol)

gene_order <- c(target_order, extra_order)
n_right <- length(gene_order)
right_rank <- setNames(seq_along(gene_order), gene_order)

nodes_left <- data.frame(
  id = mirna_order, x = 0, y = seq(1, n_right, length.out = n_left),
  label = mirna_order,
  direction = ifelse(mirna_logfc > 0, "Up", "Down"),
  stringsAsFactors = FALSE
)
left_y <- setNames(nodes_left$y, nodes_left$id)

fmt_val <- function(x) {
  r <- round(x, 1)
  ifelse(r == 0, "0.0", sprintf("%.1f", r))
}
sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

nodes_right <- data.frame(
  id = gene_order, x = 1.4, y = right_rank[gene_order],
  gene_symbol = gene_order,
  logFC = c(target_key$logFC_target[match(target_order, target_key$target_symbol)],
            top25$logFC[match(extra_order, top25$gene_symbol)]),
  adj.P.Val = c(target_key$adj.P.Val_target[match(target_order, target_key$target_symbol)],
                top25$adj.P.Val[match(extra_order, top25$gene_symbol)]),
  stringsAsFactors = FALSE
) %>%
  mutate(label = paste0(gene_symbol, "  ", fmt_val(logFC), sig_stars(adj.P.Val)))

edge_df <- edges %>%
  mutate(x = 0, xend = 1.4, y = left_y[miRNA], yend = right_rank[target_symbol])

cap <- ceiling(max(abs(nodes_right$logFC)))
logfc_colors <- c("#08306B", "#4292C6", "#F7F7F7", "#EF6548", "#99000D")

network_plot <- ggplot() +
  # --- divider between the validated-target block and the top25-only block ---
  geom_segment(aes(x = -0.3, xend = 1.8, y = n_target + 0.5, yend = n_target + 0.5),
               color = "grey75", linewidth = 0.4, linetype = "22") +

  geom_segment(data = edge_df,
               aes(x = x, xend = xend, y = y, yend = yend, color = directionality_consistent),
               linewidth = 1.3, alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = "#1E6B6B", `FALSE` = "#DD8A26"), guide = "none") +

  # --- miRNA nodes (left): categorical up/down ---
  geom_point(data = nodes_left, aes(x = x, y = y, fill = direction),
             shape = 22, size = 10, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("Up" = "#99000D", "Down" = "#08306B"), guide = "none") +
  geom_text(data = nodes_left, aes(x = x - 0.05, y = y, label = label), hjust = 1, size = 6, fontface = "bold") +

  ggnewscale::new_scale_fill() +

  # --- gene nodes (right): bubble fill = logFC; fixed size (significance is
  # already shown via the label's asterisks, so a separate size legend for
  # it is redundant) ---
  geom_point(data = nodes_right, aes(x = x, y = y, fill = logFC), size = 5.5, shape = 21, color = "white", stroke = 0.6) +
  scale_fill_gradientn(
    colors = logfc_colors, values = scales::rescale(c(-cap, -cap/2, 0, cap/2, cap)),
    limits = c(-cap, cap), oob = scales::squish, guide = "none"
  ) +
  geom_text(data = nodes_right, aes(x = x + 0.1, y = y, label = label), hjust = 0, size = 5.5, fontface = "bold") +

  scale_x_continuous(limits = c(-0.6, 1.95), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, n_right + 1), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 10, 15, 30)
  )

# Legend built with ComplexHeatmap::Legend()/packLegend(), stacked
# vertically down the right side. Font/swatch sizes here are ~2x the
# marker-panel plot's (13/11pt, 4.5mm): this canvas is 17.5in tall vs. that
# plot's 8.5in for similar content density, so matching point sizes would
# render at roughly half the visual size once both are scaled to the same
# height in a multi-panel layout.
edge_legend <- Legend(
  labels = c("Inverse (consistent with repression)", "Same direction (not direct repression)"),
  legend_gp = gpar(col = c("#1E6B6B", "#DD8A26"), lwd = 6), type = "lines",
  title = "Edge directionality", title_gp = gpar(fontsize = 26, fontface = "bold"),
  labels_gp = gpar(fontsize = 22), grid_height = unit(9, "mm"), grid_width = unit(14, "mm")
)
mirna_dir_legend <- Legend(
  labels = c("Up", "Down"), legend_gp = gpar(fill = c("#99000D", "#08306B")),
  title = "miRNA direction", title_gp = gpar(fontsize = 26, fontface = "bold"),
  labels_gp = gpar(fontsize = 22), grid_height = unit(9, "mm"), grid_width = unit(9, "mm")
)
log2fc_col_fun <- circlize::colorRamp2(
  c(-cap, -cap / 2, 0, cap / 2, cap), logfc_colors
)
log2fc_legend <- Legend(
  col_fun = log2fc_col_fun, title = "Log2FC", at = seq(-cap, cap, by = cap / 2),
  direction = "vertical", legend_height = unit(5.2, "cm"), legend_width = unit(1.1, "cm"),
  title_gp = gpar(fontsize = 26, fontface = "bold"), labels_gp = gpar(fontsize = 22)
)
legend_grob <- grid.grabExpr(
  draw(packLegend(edge_legend, mirna_dir_legend, log2fc_legend, direction = "vertical", gap = unit(16, "mm")))
)

final_plot <- cowplot::plot_grid(network_plot, legend_grob, ncol = 2, rel_widths = c(2.85, 1))

ggsave(file.path(tif_dir, "supp_figure_7e_top25_mirna_combined.pdf"), plot = final_plot,
       width = 23, height = 17.5, dpi = 600, bg = "white", device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_7e_top25_mirna_combined.pdf")), plot = final_plot,
       width = 23, height = 17.5, dpi = 600, bg = "white", device = "tiff")
message("Saved: ", file.path(tif_dir, "supp_figure_7e_top25_mirna_combined.pdf"))
