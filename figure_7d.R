# =============================================================================
# Figure 7d – DEA protein Venn diagram (vs PBS Control contrasts)
# =============================================================================
#
# Description: Venn diagram of significant proteins (FDR < 0.05) across
#   contrasts vs PBS control. Requires preprocessing steps 1–2.
#
# Protein sets depend on BioMart Entrez ID mapping (step 2), which updates
#   periodically and would shift set membership if recomputed. This figure
#   instead uses the frozen gene-set membership from the original 2024
#   analysis (see figure_7_input_data/original_2024_objects/README.md) and
#   redraws it fresh with ggvenn, since the original serialized plot object
#   no longer renders correctly under the current ggplot2 version.
#
# Prerequisites: figure_7_input_data/original_2024_objects/figure_7d_venn_gene_sets_original.csv
#
# Output: figure_7d_venn_plot_dea_pro_HA.pdf
# =============================================================================
cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(ggplot2)
  library(ggvenn)
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}
OUTPUT_DIR <- "figure_7_output_data"
ORIGINAL_DIR <- "figure_7_input_data/original_2024_objects"

venn_data <- read.csv(file.path(ORIGINAL_DIR, "figure_7d_venn_gene_sets_original.csv"), stringsAsFactors = FALSE)

# Original contrast order/labels, from old_code/COAST_Pro_HA_6_Network_Analysis.R
contrast_order <- c(
  "COVID-19 (+)  mTBI (+) ADEs vs  PBS Control",
  "ODN-D-SLO3 vs  PBS Control",
  "Poly(A:U) vs  PBS Control",
  "TNF-α vs  PBS Control"
)
contrast_labels <- c(
  "COVID-19 (+) \nmTBI (+) ADEs vs \nPBS Control",
  "ODN-D-SLO3 vs \nPBS Control",
  "Poly(A:U) vs \nPBS Control",
  "TNF-α vs \nPBS Control"
)
venn_list <- setNames(
  lapply(contrast_order, function(cn) unique(venn_data$entrez_id[venn_data$contrast == cn])),
  contrast_labels
)

venn_plot_dea_pro_HA <- ggvenn(venn_list, fill_color = c("#A2B9FC", "#E86D3D", "#DF4416", "#739D52"),
       stroke_size = 0.8, set_name_size = 12,
       stroke_color = "white", stroke_alpha = 0.8,
       text_color = "black",
       text_size = 16)
venn_plot_dea_pro_HA

dir.create(OUTPUT_DIR, showWarnings = FALSE)
venn_plot_data <- aggregate(entrez_id ~ contrast, data = venn_data, FUN = length)
names(venn_plot_data) <- c("set", "n_proteins")
write.csv(venn_plot_data, file.path(OUTPUT_DIR, "figure_7d_venn_dea_data.csv"), row.names = FALSE)

path_output <- "figure_7d_venn_plot_dea_pro_HA.pdf"
ggsave(path_output, venn_plot_dea_pro_HA, device = grDevices::cairo_pdf, width = 30, height = 25, dpi = 600)
ggsave(sub("\\.pdf$", ".tif", path_output), venn_plot_dea_pro_HA, device = "tiff", width = 30, height = 25, dpi = 600)
