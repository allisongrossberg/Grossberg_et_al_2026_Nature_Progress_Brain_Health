# =============================================================================
# Supplementary Figure 7b - Reactive astrocyte marker panel (Escartin et al.)
# =============================================================================
#
# Description: Lollipop plot of significant COV_TBI_vs_Control proteins
#   that are established reactive-astrocyte markers, ranked by logFC:
#   (1) the 10 markers from Escartin et al. 2021 Table 1 ("Potential markers
#       of reactive astrocytes", Nat Neurosci, PMID 33589835) that are
#       detected and significant here, tagged with a star;
#   (2) 5 additional literature markers not in Table 1: GBP2 (PMID
#       28099414), ICAM1 (PMID 25455510), SPARCL1/Hevin (PMID 21788491),
#       NRG1 (PMID 28456012), CLU/Clusterin (PMID 41694601).
#   A colored dot tags each marker's functional category (antiviral/immune,
#   stress & chaperone, neuronal support, homeostatic transport,
#   metabolism, identity & cytoskeleton).
#
# Prerequisites: supp_figure_7b_output_data/supp_figure_7b_marker_panel_screen_results.csv
#   (from supp_figure_7b_marker_panel_screen.R);
#   supp_figure_7b_output_data/supp_figure_7b_top25_significant_results.csv
#   (from supp_figure_7b_top25_significant_results.R)
#
# Output: supp_figure_7b_marker_panel_heatmap.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
if (basename(getwd()) == "supp_figure_7") setwd("..")

output_dir <- "supp_figure_7/supp_figure_7_output_data"
tif_dir    <- "supp_figure_7"

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggstar)
})

screen_path <- file.path(output_dir, "supp_figure_7b_marker_panel_screen_results.csv")
if (!file.exists(screen_path)) stop("Run supp_figure_7b_marker_panel_screen.R first. Expected: ", screen_path)
screened <- read.csv(screen_path, stringsAsFactors = FALSE)

# Curated literature markers (Escartin Table 1 + other literature), detected
# and significant
curated <- screened %>%
  filter(detected_in_proteomics == TRUE, !is.na(adj.P.Val_COV_TBI), adj.P.Val_COV_TBI < 0.05) %>%
  dplyr::transmute(
    gene_symbol, logFC = logFC_COV_TBI, adj.P.Val = adj.P.Val_COV_TBI,
    source = source, category = category
  )

# Avoid "-0.0" from rounding small negative values down to zero
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

# Six mechanistic categories, assigned directly by gene symbol.
category_colors <- c(
  "Antiviral/immune"        = "#C0392B",
  "Stress & chaperone"      = "#8E44AD",
  "Neuronal support"        = "#1E8449",
  "Homeostatic transport"   = "#2471A3",
  "Metabolism"              = "#CA6F1E",
  "Identity & cytoskeleton" = "#616A6B"
)
category_map <- c(
  GBP2 = "Antiviral/immune", ICAM1 = "Antiviral/immune",
  HSPB1 = "Stress & chaperone", CRYAB = "Stress & chaperone", CLU = "Stress & chaperone",
  SPARCL1 = "Neuronal support", THBS1 = "Neuronal support", NRG1 = "Neuronal support",
  SLC1A2 = "Homeostatic transport", SLC1A3 = "Homeostatic transport", MT2A = "Homeostatic transport",
  FABP7 = "Metabolism", MAOB = "Metabolism",
  SOX9 = "Identity & cytoskeleton", VIM = "Identity & cytoskeleton"
)

plot_df <- curated %>%
  mutate(
    category = unname(category_map[gene_symbol]),
    dot_color = unname(ifelse(!is.na(category), category_colors[category], NA)),
    has_star = source == "Escartin Table 1",
    label = paste0(fmt_val(logFC), sig_stars(adj.P.Val)),
    label_hjust = ifelse(logFC >= 0, -0.35, 1.35)
  )

# Order top (up) to bottom (down) by logFC. gene_symbol (plain text, no
# embedded symbols) is the shared y factor for both the lollipop plot and
# the symbol column below, so their rows line up exactly.
y_order <- plot_df %>% arrange(desc(logFC)) %>% pull(gene_symbol)
plot_df <- plot_df %>% mutate(gene_symbol = factor(gene_symbol, levels = rev(y_order)))

cap <- ceiling(max(abs(plot_df$logFC), na.rm = TRUE))

lollipop_plot <- ggplot(plot_df, aes(x = logFC, y = gene_symbol)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
  geom_segment(aes(xend = 0, yend = gene_symbol, color = logFC), linewidth = 2.0) +
  geom_point(aes(color = logFC), size = 4) +
  geom_text(aes(label = label, hjust = label_hjust), color = "black", fontface = "bold", size = 4.4) +
  scale_color_gradientn(
    colors = c("#08306B", "#4292C6", "#F7F7F7", "#EF6548", "#99000D"),
    values = scales::rescale(c(-cap, -cap/2, 0, cap/2, cap)),
    limits = c(-cap, cap), guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.32, 0.32))) +
  labs(x = "Log2 Fold Change: COVID-19+mTBI+ vs COVID-19-mTBI-", y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 13, face = "bold", color = "gray15"),
    axis.text.x = element_text(size = 13, face = "bold", color = "gray15"),
    axis.title.x = element_text(size = 14, face = "bold", color = "gray15", margin = margin(t = 12, b = 6)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 10, 6, 8)
  )

# Symbols as a narrow aligned column to the left of the gene names, using
# the same y factor as the main plot and cowplot::plot_grid(align = "h") so
# the rows line up exactly. Every layer here uses the full plot_df (all 15
# rows), not a row-subset: two layers with different row subsets of the
# same shared factor would each build their discrete-position mapping from
# only the levels present in that layer, scrambling y-positions for genes
# missing from one layer but not the other. Absence is encoded via a
# transparent color / blank label instead.
symbol_plot <- ggplot(plot_df, aes(y = gene_symbol)) +
  geom_point(aes(x = 0, color = I(ifelse(is.na(dot_color), "transparent", dot_color))),
             size = 3.2) +
  # ggstar::geom_star draws a vector polygon rather than a Unicode glyph,
  # which avoids font/encoding issues when this PDF is placed into Adobe
  # Illustrator without linking (see full note below).
  ggstar::geom_star(aes(x = 1, fill = I(ifelse(has_star, "#B8860B", "transparent"))),
                     starshape = 1, colour = NA, size = 3.2) +
  scale_x_continuous(limits = c(-0.6, 1.6), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(10, 0, 10, 2))

# Legend built with ComplexHeatmap::Legend()/packLegend(), stacked
# vertically down the right side.
log2fc_col_fun <- circlize::colorRamp2(
  c(-cap, -cap / 2, 0, cap / 2, cap),
  c("#08306B", "#4292C6", "#F7F7F7", "#EF6548", "#99000D")
)
log2fc_legend <- Legend(
  col_fun = log2fc_col_fun, title = "Log2FC", at = c(-cap, -cap / 2, 0, cap / 2, cap),
  direction = "vertical", legend_height = unit(2.6, "cm"),
  title_gp = gpar(fontsize = 13, fontface = "bold"), labels_gp = gpar(fontsize = 11)
)
category_legend <- Legend(
  labels = names(category_colors), legend_gp = gpar(fill = unname(category_colors)),
  title = "Category", title_gp = gpar(fontsize = 13, fontface = "bold"),
  labels_gp = gpar(fontsize = 11), grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")
)
# Legend key uses the same vector-drawn star as the plot's ggstar::geom_star
# markers (see note above), not a Unicode "★" glyph or pch = 8.
draw_star_key <- function(x, y, w, h) {
  sz <- min(as.numeric(w), as.numeric(h)) * 0.8
  grid.draw(ggstar:::starGrob(x = x, y = y, starshape = 1, angle = 0,
                               gp = gpar(fill = "#B8860B", col = NA, fontsize = sz),
                               position.units = "mm", size.units = "mm"))
}
star_legend <- Legend(
  labels = "Escartin et al. 2021\nNat Neurosci Table 1\n(PMID 33589835)",
  type = "grid", graphics = list(draw_star_key),
  size = unit(4, "mm"), labels_gp = gpar(fontsize = 10.5), title = "Key",
  title_gp = gpar(fontsize = 13, fontface = "bold")
)
legend_grob <- grid.grabExpr(
  draw(packLegend(log2fc_legend, category_legend, star_legend, direction = "vertical", gap = unit(8, "mm")))
)

main_row <- cowplot::plot_grid(
  symbol_plot, lollipop_plot,
  ncol = 2, rel_widths = c(0.45, 4), align = "h", axis = "tb"
)
final_plot <- cowplot::plot_grid(main_row, legend_grob, ncol = 2, rel_widths = c(4.3, 1.5))

dir.create(tif_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(tif_dir, "supp_figure_7b_marker_panel_heatmap.pdf"), plot = final_plot,
       width = 8.1, height = 8.5, dpi = 600, bg = "white", device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_7b_marker_panel_heatmap.pdf")), plot = final_plot,
       width = 8.1, height = 8.5, dpi = 600, bg = "white", device = "tiff")
message("Saved: ", file.path(tif_dir, "supp_figure_7b_marker_panel_heatmap.pdf"))
