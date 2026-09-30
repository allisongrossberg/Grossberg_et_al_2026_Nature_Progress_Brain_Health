# =============================================================================
# Figure 7c – Circos/chord diagram: GSEA KEGG pathways by contrast
# =============================================================================
#
# Description: Chord diagram linking contrasts to shared KEGG pathways (from
#   GSEA). Requires preprocessing steps 1–4 and pro_HA_GSEA_KEGG_*.csv files.
#
# This plot depends on GSEA KEGG results (step 4), which query the KEGG
#   REST API; pathway names and membership change between KEGG releases, so
#   recomputing from current KEGG results would not reproduce the published
#   diagram. It instead loads the frozen condition x pathway enrichment
#   matrix from the original 2024 analysis workspace (see
#   figure_7_input_data/original_2024_objects/README.md for provenance).
#
# Prerequisites: figure_7_input_data/original_2024_objects/circos_matrix.rds
#
# Output: figure_7c_circos_plot.pdf
# =============================================================================
cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
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

library(circlize)
library(dplyr)
library(stringr)

circos_matrix <- as.matrix(readRDS(file.path(ORIGINAL_DIR, "circos_matrix.rds")))

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Define a named vector for condition renaming
condition_rename <- c(
  "COV_TBI_vs_NegCTRL" = "COVID-19 (+) mTBI (+) ADEs vs PBS Control",
  "COV_TBI_vs_Control" = "COVID-19 (+) mTBI (+) ADEs vs COVID-19 (-) mTBI (-) ADEs",
  "Poly_A_U_vs_NegCTRL" = "Poly A:U vs PBS Control",
  "Poly_A_U_vs_Control" = "Poly A:U vs COVID-19 (-) mTBI (-) ADEs",
  "TNFa_vs_NegCTRL" = "TNF-α vs PBS Control",
  "TNFa_vs_Control" = "TNF-α vs COVID-19 (-) mTBI (-) ADEs",
  "IL1b_vs_NegCTRL" = "IL-1β vs PBS Control",
  "IL1b_vs_Control" = "IL-1β vs COVID-19 (-) mTBI (-) ADEs",
  "LPS_vs_NegCTRL" = "LPS vs PBS Control",
  "LPS_vs_Control" = "LPS vs COVID-19 (-) mTBI (-) ADEs",
  "ODN_D_SLO3_vs_NegCTRL" = "ODN-D-SLO3 vs COVID-19 (-) mTBI (-) ADEs",
  "ODN_D_SLO3_vs_Control" = "ODN-D-SLO3 vs PBS Control"
)

conditions <- rownames(circos_matrix)
pathways <- colnames(circos_matrix)

# Create a function to generate more distinct colors
get_distinct_colors <- function(n) {
  hues = seq(15, 375, length = n + 1)
  colors = hcl(h = hues, l = 65, c = 100)[1:n]
  colors = colorRampPalette(colors)(n)
  return(colors)
}

# Generate colors
condition_colors <- get_distinct_colors(length(conditions))
pathway_colors <- get_distinct_colors(length(pathways))

# Create the named vector for grid.col
grid_colors <- c(condition_colors, pathway_colors)
names(grid_colors) <- c(conditions, pathways)

# Calculate text size with a larger base size
text_size <- 0.8 * min(1, 15 / max(nchar(c(condition_rename, pathways))))
text_size <- max(text_size, 0.99)  # Ensure the text size is not too small

add_line_breaks <- function(label, max_width = 20) {
  words <- strsplit(label, " ")[[1]]
  lines <- character(0)
  current_line <- words[1]

  for (word in words[-1]) {
    if (nchar(current_line) + nchar(word) + 1 > max_width) {
      lines <- c(lines, current_line)
      current_line <- word
    } else {
      current_line <- paste(current_line, word)
    }
  }
  lines <- c(lines, current_line)

  paste(lines, collapse = "\n")
}

condition_rename <- sapply(condition_rename, add_line_breaks)

draw_circos_plot <- function() {
  # Adjust circos parameters
  circos.par(gap.after = c(rep(2, length(conditions) - 1), 8, rep(1, length(pathways) - 1), 8),
             start.degree = 90,
             track.height = 0.15)

  # Create the Circos plot
  par(mar=c(1,1,1,1))
  circos.clear()

  chordDiagram(circos_matrix,
               grid.col = grid_colors,
               transparency = 0.5,
               annotationTrack = "grid",
               preAllocateTracks = 1)

  # Add labels with improved visibility and larger text
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.index <- CELL_META$sector.index
    label <- condition_rename[sector.index]
    if (is.na(label)) label <- add_line_breaks(sector.index)
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], label,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = text_size, col = "black", font = 2)
  }, bg.border = NA)
}

# Canvas: 6500x4000px @ 300dpi = 21.67x13.33in; circos panel gets 4/7 of the width.
path_output <- "figure_7c_circos_plot.pdf"
circos_width <- 21.666667 * 4 / 7
circos_height <- 13.333333
grDevices::cairo_pdf(path_output, width = circos_width, height = circos_height)
draw_circos_plot()
grDevices::dev.off()
message("Saved: ", path_output)

tif_output <- sub("\\.pdf$", ".tif", path_output)
grDevices::tiff(tif_output, width = circos_width, height = circos_height, units = "in", res = 600, compression = "lzw")
draw_circos_plot()
grDevices::dev.off()
message("Saved: ", tif_output)
