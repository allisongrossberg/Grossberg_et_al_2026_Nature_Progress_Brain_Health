# =============================================================================
# Figure 7i – Protein interaction network: Poly A:U vs PBS Control
# =============================================================================
#
# Description: Network plot for Poly A:U vs Neg control. Requires steps 1–5.
#
# Node/edge data depends on GSEA KEGG results (step 4), which query the
#   live KEGG API, so re-running from scratch would not exactly reproduce
#   this network. This figure instead uses the frozen igraph object from
#   the original 2024 analysis (see figure_7_input_data/
#   original_2024_objects/README.md) and redraws it fresh with ggraph,
#   since the original serialized plot object no longer renders correctly
#   under the current ggplot2 version.
#
# Prerequisites: figure_7_input_data/original_2024_objects/poly_a_u_neg_control_plot.rds
#
# Output: figure_7i_Poly_AU_vs_Neg_Cont_network_plot.pdf
# =============================================================================
cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
OUTPUT_DIR <- "figure_7_output_data"
ORIGINAL_DIR <- "figure_7_input_data/original_2024_objects"
library(ggplot2)
library(igraph)
library(ggraph)
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}

poly_a_u_neg_control_plot <- readRDS(file.path(ORIGINAL_DIR, "poly_a_u_neg_control_plot.rds"))

dir.create(OUTPUT_DIR, showWarnings = FALSE)
g <- poly_a_u_neg_control_plot$g
ed <- as.data.frame(igraph::as_edgelist(g, names = TRUE))
names(ed) <- c("V1", "V2")
write.csv(ed, file.path(OUTPUT_DIR, "figure_7i_network_edges.csv"), row.names = FALSE)
write.csv(ed, file.path(OUTPUT_DIR, "pro_HA_network_edges_Poly_A_U_vs_NegCTRL.csv"), row.names = FALSE)
if (igraph::vcount(g) > 0) {
  write.csv(data.frame(node = igraph::V(g)$name, gene_name = igraph::V(g)$gene_name, pathway = igraph::V(g)$pathway, direction = igraph::V(g)$direction, degree = igraph::degree(g)), file.path(OUTPUT_DIR, "figure_7i_network_nodes.csv"), row.names = FALSE)
}

# Rebuilt fresh from the frozen igraph object; same call as the original
# generate_proteomics_network_plot() (old_code/COAST_Pro_HA_6_Network_Analysis.R).
custom_colors <- c("#1E3221", "#395724", "#739D52", "#A2B9FC",
                   "#7583B7", "#DF4416",  "#E86D3D", "#FFD35A",
                   "#FCDC94", "#FEFBD8", "#EECEB9",  "#BB9AB1",
                   "#D66EB7", "#8C3061", "#522258", "#201E43")
build_network_plot <- function(g) {
  set.seed(123)
  ggraph(g, layout = "stress") +
    geom_edge_link(aes(edge_alpha = combined_score, edge_width = combined_score),
                   edge_colour = "#B4B4B3") +
    geom_node_point(aes(size = degree, color = pathway, shape = node_shape)) +
    geom_node_text(aes(label = gene_name), repel = TRUE, size = 4.5,
                   max.overlaps = 15) +
    scale_edge_width_continuous(range = c(0.1, 0.7), name = "Interaction Strength") +
    scale_size_continuous("Degree", range = c(2, 8)) +
    scale_color_manual(values = custom_colors, name = "Pathway") +
    theme_graph(base_family = "Helvetica") +
    theme(legend.position = "right",
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22)) +
    guides(size = guide_legend(order = 1, title = "Node Degree"),
           edge_width = guide_legend(order = 2, title = "Interaction Strength"),
           edge_alpha = "none",
           color = guide_legend(order = 3, title = "Pathway",
                                override.aes = list(size = 5)),
           shape = guide_legend(order = 4, title = "Regulation",
                                override.aes = list(size = 5))
           )
}

path_out <- file.path(getwd(), "figure_7i_Poly_AU_vs_Neg_Cont_network_plot.pdf")
network_plot <- build_network_plot(g)

grDevices::cairo_pdf(path_out, width = 20, height = 12)
print(network_plot)
grDevices::dev.off()
message("Saved: ", path_out)

tif_out <- sub("\\.pdf$", ".tif", path_out)
grDevices::tiff(tif_out, width = 20, height = 12, units = "in", res = 600, compression = "lzw")
print(network_plot)
grDevices::dev.off()
message("Saved: ", tif_out)
