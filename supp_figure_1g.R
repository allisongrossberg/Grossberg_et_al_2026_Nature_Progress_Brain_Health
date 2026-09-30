# =============================================================================
# Supplementary Figure 1g – Biomarkers by time since exposure (COVID-19 vs. mTBI)
# =============================================================================
#
# Description: Combines the two "time since exposure" analyses
#   (supp_figure_1g_preprocessing_covid_recency.R,
#   supp_figure_1g_preprocessing_tbi_recency.R) into a single 2x2 figure:
#   biomarker (row) x exposure type (column), showing the pooled
#   ("Overall") Kruskal-Wallis comparison across Recent/Intermediate/Remote
#   tertiles for each. Per-subgroup breakdowns and Dunn post-hoc
#   comparisons remain in the individual preprocessing scripts' own
#   figures and are not duplicated here.
#
# Prerequisites: supp_figure_1g_preprocessing_covid_recency.R and
#   supp_figure_1g_preprocessing_tbi_recency.R must be run first (this
#   script reads their per-participant output CSVs).
#
# Output: supp_figure_1g_panel_recency.pdf; supp_figure_1g_panel_recency_stats.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
# Normalize working directory to repo root regardless of where this is run from.
if (basename(getwd()) == "supp_figure_1") setwd("..")
output_dir <- "supp_figure_1/supp_figure_1_output_data"
tif_dir    <- "supp_figure_1"

path_covid <- file.path(output_dir, "supp_figure_1_covid_recency_input_data.csv")
path_tbi   <- file.path(output_dir, "supp_figure_1_tbi_recency_input_data.csv")
missing <- c(
  if (!file.exists(path_covid)) "supp_figure_1_covid_recency_input_data.csv (run supp_figure_1g_preprocessing_covid_recency.R)",
  if (!file.exists(path_tbi))   "supp_figure_1_tbi_recency_input_data.csv (run supp_figure_1g_preprocessing_tbi_recency.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

covid_df <- read.csv(path_covid, stringsAsFactors = FALSE) %>%
  dplyr::mutate(exposure = "COVID-19")
tbi_df <- read.csv(path_tbi, stringsAsFactors = FALSE) %>%
  dplyr::mutate(exposure = "mTBI")

Combined_Data <- dplyr::bind_rows(covid_df, tbi_df)
Combined_Data$recency_category <- factor(Combined_Data$recency_category,
                                          levels = c("Recent", "Intermediate", "Remote"))
Combined_Data$exposure <- factor(Combined_Data$exposure, levels = c("COVID-19", "mTBI"))

# Long format: one row per participant x biomarker
plot_data <- Combined_Data %>%
  tidyr::pivot_longer(cols = c(Abeta42_Abeta40_Ratio_R, vib_anti_abeta_1_42_igg_iga),
                       names_to = "biomarker_name", values_to = "biomarker_value") %>%
  dplyr::mutate(biomarker_name = dplyr::recode(biomarker_name,
                                                Abeta42_Abeta40_Ratio_R = "Aβ42/Aβ40",
                                                vib_anti_abeta_1_42_igg_iga = "Anti-Aβ IgG/IgA")) %>%
  dplyr::filter(!is.na(biomarker_value), !is.na(recency_category))

# Kruskal-Wallis omnibus per exposure x biomarker panel (pooled/Overall), FDR across the 4 panels
panels <- plot_data %>% dplyr::distinct(exposure, biomarker_name)
panel_stats <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  ex <- panels$exposure[i]; b <- panels$biomarker_name[i]
  sub <- plot_data %>% dplyr::filter(exposure == ex, biomarker_name == b)
  test <- kruskal.test(biomarker_value ~ droplevels(recency_category), data = sub)
  data.frame(exposure = ex, biomarker_name = b, n = nrow(sub), kruskal_p = test$p.value)
}))
panel_stats$p_fdr <- p.adjust(panel_stats$kruskal_p, method = "fdr")
write.csv(panel_stats, file.path(output_dir, "supp_figure_1g_panel_recency_stats.csv"), row.names = FALSE)
cat("Panel B (recency) stats:\n"); print(panel_stats)

# The "COVID-19" and "mTBI" facets are overlapping stratifications, not
# separate populations: the double-exposed group appears in both (see
# supp_figure_1f.R for the exact overlap count). Points are colored by
# qq_group (matching figure_1k.R's radar-plot color mapping) to make that
# overlap visible rather than coloring by recency_category, which is
# redundant with the x-axis position. Only 3 of the 4 qq_group levels
# appear here, since neither facet includes the double-negative group.
qq_group_colors <- c(
  "COVID-19 (+) mTBI (-)" = "#739D51",
  "COVID-19 (-) mTBI (+)" = "#D06FAC",
  "COVID-19 (+) mTBI (+)" = "#DD4726"
)

# One row per biomarker (stacked via cowplot) so each row can carry its own
# y-axis label -- a shared facet_grid can only show a single y-axis title
# for both rows. Box/jitter parameters, palette, and theme match
# figure_1i.R / figure_1j.R.
make_biomarker_row <- function(biomarker, y_label, show_x, show_strip, show_legend) {
  ggplot(plot_data %>% dplyr::filter(biomarker_name == biomarker),
         aes(x = recency_category, y = biomarker_value)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.5, fill = "white", color = "black",
                 alpha = 0.7, outlier.shape = NA, linewidth = 0.8) +
    geom_jitter(aes(fill = qq_group), width = 0.2, height = 0, shape = 21, size = 4, stroke = 0.5) +
    scale_fill_manual(values = qq_group_colors, name = NULL) +
    facet_wrap(~ exposure, nrow = 1, scales = "free_y") +
    theme_minimal() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(size = 14, family = "Arial", face = "bold", colour = "black"),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      strip.text = if (show_strip) element_text(size = 16, lineheight = .9, family = "Arial", face = "bold", colour = "black") else element_blank(),
      axis.text.x = if (show_x) element_text(size = 16, lineheight = .9, family = "Arial", face = "bold", colour = "black") else element_blank(),
      axis.title.x = if (show_x) element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)) else element_blank(),
      axis.text.y = element_text(size = 16, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
      axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0))
    ) +
    labs(x = if (show_x) "Time Since Incidence (Tertile)" else NULL, y = y_label)
}

# Legend only on the bottom row, so the combined figure gets one shared
# legend instead of two.
igg_row   <- make_biomarker_row("Anti-Aβ IgG/IgA", "Anti-Aβ1-42 (IgG/IgA)", show_x = FALSE, show_strip = TRUE, show_legend = FALSE)
abeta_row <- make_biomarker_row("Aβ42/Aβ40", "Plasma Aβ42/Aβ40 Ratio", show_x = TRUE, show_strip = FALSE, show_legend = TRUE)

Panel_Recency_Plot <- cowplot::plot_grid(igg_row, abeta_row, ncol = 1, align = "v", axis = "lr")

ggsave(file.path(tif_dir, "supp_figure_1g_panel_recency.pdf"), Panel_Recency_Plot, width = 11, height = 10, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1g_panel_recency.pdf")), Panel_Recency_Plot, width = 11, height = 10, dpi = 600, device = "tiff")
