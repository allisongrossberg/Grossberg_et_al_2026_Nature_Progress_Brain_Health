# =============================================================================
# Figure 1k – Radar chart of summary scores by group
# =============================================================================
#
# Description: Radar chart of normalized outcome scores (Neuro-QoL, N3PA,
#   anti-Aβ IgG/IgA, etc.) by group. Run preprocessing and figure_1f, 1i, 1j first.
#
#   Every axis is recoded before scaling so that outward consistently means
#   worse outcome / greater burden, regardless of the raw score's clinical
#   direction.
#
# Prerequisites: figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv,
#   figure_1f_NeuroQOL_plot_data.csv, figure_1i_N3PA_Abeta_plot_data.csv,
#   figure_1j_Vib_Figures_Data.csv (run preprocessing then figure_1f, 1i, 1j first).
#
# Output: figure_1k_radar_plot.pdf; data in figure_1_output_data
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
output_dir <- if (dir.exists("figure_1_revised")) "figure_1_revised/figure_1_output_data" else "figure_1_output_data"
tif_dir    <- if (dir.exists("figure_1_revised")) "figure_1_revised" else "."
library(dplyr)
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(
  file.path(output_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"),
  stringsAsFactors = FALSE
)
NeuroQOL_Data_join <- read.csv(file.path(output_dir, "figure_1f_NeuroQOL_plot_data.csv"), stringsAsFactors = FALSE)
N3PA_plot_data    <- read.csv(file.path(output_dir, "figure_1i_N3PA_Abeta_plot_data.csv"), stringsAsFactors = FALSE)
Vib_Figures_Data  <- read.csv(file.path(output_dir, "figure_1j_Vib_Figures_Data.csv"), stringsAsFactors = FALSE)

COAST_Study_Data_Combined_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  select(participant_id, qq_group, age_years, qq_phq8_average_score, qq_eq5d_index_score, qq_wai_2,
         recode_overall_neuro_psych_sev_score, recode_overall_neuro_psych_freq_score) %>%
  left_join(NeuroQOL_Data_join %>% select(participant_id, qq_group, Neuro_QOL_TScore),
            by = c("participant_id", "qq_group")) %>%
  left_join(N3PA_plot_data %>% select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R),
            by = c("participant_id", "qq_group")) %>%
  left_join(Vib_Figures_Data %>% select(participant_id, qq_group, vib_anti_abeta_1_42_igg_iga),
            by = c("participant_id", "qq_group"))

radar_cols <- c("participant_id", "qq_group", "qq_phq8_average_score", "qq_eq5d_index_score", "qq_wai_2",
                "recode_overall_neuro_psych_sev_score", "recode_overall_neuro_psych_freq_score",
                "Neuro_QOL_TScore", "Abeta42_Abeta40_Ratio_R", "vib_anti_abeta_1_42_igg_iga")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(COAST_Study_Data_Combined_Data %>% dplyr::select(dplyr::all_of(radar_cols[radar_cols %in% names(COAST_Study_Data_Combined_Data)])), file.path(output_dir, "figure_1k_radar_plot_data.csv"), row.names = FALSE)

library(fmsb)
library(dplyr)
library(tibble)

# Group means. Columns are ordered into three adjacent clusters (symptom/
# psychiatric burden, functional/quality-of-life, biomarkers) so related
# axes sit next to each other around the wheel.
radar_means <- COAST_Study_Data_Combined_Data %>%
  dplyr::group_by(qq_group) %>%
  dplyr::summarise(across(c(qq_phq8_average_score, recode_overall_neuro_psych_sev_score, recode_overall_neuro_psych_freq_score,
                     qq_eq5d_index_score, qq_wai_2, Neuro_QOL_TScore,
                     Abeta42_Abeta40_Ratio_R, vib_anti_abeta_1_42_igg_iga),
                   mean, na.rm = TRUE)) %>%
  column_to_rownames("qq_group")

# Column labels (same order as above).
clean_names <- c(
  "PHQ-8 Score",
  "Neuro/Psych\nSeverity",
  "Neuro/Psych\nFrequency",
  "EQ-5D Score",
  "WAI Score",
  "Neuro-QOL",
  "Aβ42/Aβ40",
  "Anti-Aβ IgG/IgA"
)
colnames(radar_means) <- clean_names

# Axes where a higher raw score is clinically better are flipped after
# scaling so that higher (scaled) consistently means worse.
higher_is_better <- c("EQ-5D Score", "WAI Score", "Neuro-QOL", "Aβ42/Aβ40")

# Each axis is scaled against the pooled 5th-95th percentile of all
# individual participants (across groups), not the range spanned by the
# 4 group means -- scaling to the group-mean range would force the most
# extreme group mean to the plot boundary on every axis regardless of how
# small that difference is relative to real participant-level variation.
raw_cols <- c("qq_phq8_average_score", "recode_overall_neuro_psych_sev_score", "recode_overall_neuro_psych_freq_score",
              "qq_eq5d_index_score", "qq_wai_2", "Neuro_QOL_TScore",
              "Abeta42_Abeta40_Ratio_R", "vib_anti_abeta_1_42_igg_iga")
ref_range <- sapply(COAST_Study_Data_Combined_Data[raw_cols], function(x) quantile(x, c(0.05, 0.95), na.rm = TRUE))
colnames(ref_range) <- clean_names
write.csv(data.frame(bound = rownames(ref_range), ref_range, check.names = FALSE),
          file.path(output_dir, "figure_1k_radar_plot_reference_range_p5_p95.csv"), row.names = FALSE)

radar_scaled <- as.data.frame(mapply(function(x, name) {
  lo <- ref_range["5%", name]; hi <- ref_range["95%", name]
  pmin(pmax((x - lo) / (hi - lo), 0), 1)  # clip in case a group mean falls outside [p5, p95]
}, radar_means, colnames(radar_means), SIMPLIFY = FALSE), check.names = FALSE)
rownames(radar_scaled) <- rownames(radar_means)
radar_scaled[higher_is_better] <- 1 - radar_scaled[higher_is_better]

write.csv(tibble::rownames_to_column(radar_scaled, "qq_group"),
          file.path(output_dir, "figure_1k_radar_plot_scaled_group_means.csv"), row.names = FALSE)

# Add explicit max/min rows (already scaled to 0-1, so max = 1 and min = 0
# for every axis, keeping the "outward = worse" convention consistent)
radar_data <- rbind(
  rep(1, ncol(radar_scaled)),
  rep(0, ncol(radar_scaled)),
  radar_scaled
)

# Define colors for each group
color_mapping <- c(
  "COVID-19 (+) mTBI (+)" = "#DD4726",
  "COVID-19 (-) mTBI (-)" = "#7583B7",
  "COVID-19 (+) mTBI (-)" = "#739D51",
  "COVID-19 (-) mTBI (+)" = "#D06FAC"
)

group_names <- rownames(radar_data)[-(1:2)]  # first two rows are max/min, not groups

colors <- color_mapping[group_names]

# order
desired_order <- c("COVID-19 (-) mTBI (-)", "COVID-19 (+) mTBI (-)",
                   "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)")
legend_order <- match(desired_order, names(color_mapping))
ordered_group_names <- names(color_mapping)[legend_order]
ordered_colors <- color_mapping[ordered_group_names]

draw_radar_plot <- function() {
  # Create the radar plot (larger text for readability)
  # Bottom margin holds the 2-row legend below the chart; left/right margin
  # widened so the two horizontal vertex labels (pushed out further below)
  # aren't clipped (radarchart uses ~ -1.2 to 1.2 in y)
  par(mar = c(11, 3, 1, 3), font = 2, cex.axis = 1.4)

  # Axis labels are placed by fmsb at a fixed radius (1.2) around the wheel.
  # "Neuro/Psych Frequency" (left, horizontal) and "Aβ42/Aβ40" (right,
  # horizontal) sit close enough to that radius to overlap the filled polygons,
  # so those two are blanked out here and redrawn further out below; the rest
  # keep fmsb's default placement.
  vlabels <- c("PHQ-8 Score", "Neuro/Psych\nSeverity", "", "EQ-5D Score",
               "WAI Score", "Neuro-QOL", "", "Anti-Aβ IgG/IgA")

  radarchart(
    radar_data,
    pcol = colors,
    pfcol = scales::alpha(colors, 0.5),
    plwd = 3,
    plty = 1,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey",
    caxislabels = c("Best", "", "", "", "Worst"),
    cglwd = 0.8,
    vlabels = vlabels,
    vlcex = 1.6
  )

  par(xpd = NA)

  # Manually redraw the two overlapping vertex labels a bit further from center
  # than fmsb's default 1.2 radius (angles computed the same way fmsb does)
  n_axes <- ncol(radar_data)
  theta <- seq(90, 450, length.out = n_axes + 1) * pi / 180
  theta <- theta[1:n_axes]
  label_radius <- 1.32
  text(cos(theta[3]) * label_radius, sin(theta[3]) * label_radius, "Neuro/Psych\nFrequency", cex = 1.6, font = 2)
  text(cos(theta[7]) * label_radius, sin(theta[7]) * label_radius, "Aβ42/Aβ40", cex = 1.6, font = 2)

  par(xpd = TRUE)

  # Legend below chart: 2x2 layout, centered underneath
  legend(
    x = 0,
    y = -1.45,
    legend = ordered_group_names,
    col = ordered_colors,
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 1.5,
    ncol = 2,
    horiz = FALSE,
    xjust = 0.5
  )
  par(xpd = FALSE)
}

# Wider figure so horizontal legend fits without truncation
pdf_path <- file.path(tif_dir, "figure_1k_radar_plot.pdf")
cairo_pdf(pdf_path, width = 14, height = 10)
draw_radar_plot()
dev.off()

tif_path <- sub("\\.pdf$", ".tif", pdf_path)
tiff(tif_path, width = 14, height = 10, units = "in", res = 600, compression = "lzw")
draw_radar_plot()
dev.off()
