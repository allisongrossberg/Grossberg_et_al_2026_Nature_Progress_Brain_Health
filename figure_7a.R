# =============================================================================
# Figure 7a – Number of significant proteins per contrast (faceted bar graph)
# =============================================================================
#
# Description: Faceted bar plot of counts of significant proteins (FDR < 0.05,
#   |log2FC| > 2) for each DEA contrast. Requires preprocessing steps 1 and 2.
#
# Prerequisites: figure_7_output_data (run preprocessing scripts first)
#
# Output: figure_7a_pro_HA_bar_plot.pdf
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
sig_counts_file <- file.path(OUTPUT_DIR, "pro_HA_significant_protein_counts.csv")
if (!file.exists(sig_counts_file)) {
  stop("Run figure_7_preprocessing_1_data_cleaning.R and figure_7_preprocessing_step_2_DEA.R first. Expected: ", sig_counts_file)
}
sig_counts <- read.csv(sig_counts_file, stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(stringr)

# Build summary from saved counts (same structure as before)
summary_df_pro <- data.frame(
  Contrast = sig_counts$contrast,
  Number_of_Significant_proteins = sig_counts$n_significant,
  stringsAsFactors = FALSE)

# Define the renaming vector
condition_rename <- c(
  "Cont_Exo_vs_NegCTRL" = "COVID-19 (-) mTBI (-) ADEs vs PBS Control",
  "Cont_Exo_vs_COV_TBI" = "COVID-19 (-) mTBI (-) ADEs vs COVID-19 (+) mTBI (+) ADEs",
  "COV_Only_vs_NegCTRL" = "COVID-19 (+) mTBI (-) ADEs vs PBS Control",
  "COV_Only_vs_Control" = "COVID-19 (+) mTBI (-) ADEs vs COVID-19 (-) mTBI (-) ADEs",
  "COV_Only_vs_COV_TBI" = "COVID-19 (+) mTBI (-) ADEs vs COVID-19 (+) mTBI (+) ADEs",
  "TBI_Only_vs_NegCTRL" = "COVID-19 (-) mTBI (+) ADEs vs PBS Control",
  "TBI_Only_vs_Control" = "COVID-19 (-) mTBI (+) ADEs vs COVID-19 (-) mTBI (-) ADEs",
  "TBI_Only_vs_COV_TBI" = "COVID-19 (-) mTBI (+) ADEs vs COVID-19 (+) mTBI (+) ADEs",
  "TBI_Only_vs_COV_Only" = "COVID-19 (-) mTBI (+) ADEs vs COVID-19 (+) mTBI (-) ADEs",
  "COV_TBI_vs_NegCTRL" = "COVID-19 (+) mTBI (+) ADEs vs PBS Control",
  "COV_TBI_vs_Control" = "COVID-19 (+) mTBI (+) ADEs vs COVID-19 (-) mTBI (-) ADEs",
  "Poly_A_U_vs_NegCTRL" = "Poly A:U vs PBS Control",
  "Poly_A_U_vs_Control" = "Poly A:U vs COVID-19 (-) mTBI (-) ADEs",
  "Poly_A_U_vs_COV_TBI" = "Poly A:U vs COVID-19 (+) mTBI (+) ADEs",
  "TNFa_vs_NegCTRL" = "TNF-α vs PBS Control",
  "TNFa_vs_Control" = "TNF-α vs COVID-19 (-) mTBI (-) ADEs",
  "TNFa_vs_COV_TBI" = "TNF-α vs COVID-19 (+) mTBI (+) ADEs",
  "IL1b_vs_NegCTRL" = "IL-1β vs PBS Control",
  "IL1b_vs_Control" = "IL-1β vs COVID-19 (-) mTBI (-) ADEs",
  "IL1b_vs_COV_TBI" = "IL-1β vs COVID-19 (+) mTBI (+) ADEs",
  "LPS_vs_NegCTRL" = "LPS vs PBS Control",
  "LPS_vs_Control" = "LPS vs COVID-19 (-) mTBI (-) ADEs",
  "LPS_vs_COV_TBI" = "LPS vs COVID-19 (+) mTBI (+) ADEs",
  "ODN_D_SLO3_vs_NegCTRL" = "ODN-D-SLO3 vs PBS Control",
  "ODN_D_SLO3_vs_Control" = "ODN-D-SLO3 vs COVID-19 (-) mTBI (-) ADEs",
  "ODN_D_SLO3_vs_COV_TBI" = "ODN-D-SLO3 vs COVID-19 (+) mTBI (+) ADEs"
)

# Function to categorize the contrasts
categorize_contrast <- function(contrast) {
  if (grepl("vs PBS Control$", contrast)) {
    return("vs PBS Control")
  } else if (grepl("vs COVID-19 \\(-\\) mTBI \\(-\\) ADEs$", contrast)) {
    return("vs COVID-19 (-) mTBI (-)")
  } else if (grepl("vs COVID-19 \\(\\+\\) mTBI \\(\\+\\) ADEs$", contrast)) {
    return("vs COVID-19 (+) mTBI (+)")
  } else if (grepl("vs COVID-19 \\(\\+\\) mTBI \\(-\\) ADEs$", contrast)) {
    return("vs COVID-19 (+) mTBI (-)")
  } else {
    return("Other")
  }
}

# Function to trim the contrast labels
trim_contrast_label <- function(contrast) {
  return(str_remove(contrast, " vs .*$"))
}

# Process the dataframe
summary_df_pro <- summary_df_pro %>%
  mutate(
    Contrast = case_when(
      Contrast %in% names(condition_rename) ~ condition_rename[Contrast],
      TRUE ~ Contrast
    ),
    Category = sapply(Contrast, categorize_contrast),
    Trimmed_Contrast = sapply(Contrast, trim_contrast_label)
  )

# Set the order of categories
category_order <- c("vs PBS Control", "vs COVID-19 (-) mTBI (-)", "vs COVID-19 (+) mTBI (-)", "vs COVID-19 (+) mTBI (+)")
summary_df_pro$Category <- factor(summary_df_pro$Category, levels = category_order)

# Create a unique identifier for each bar across all facets
summary_df_pro <- summary_df_pro %>%
  arrange(Category, Trimmed_Contrast) %>%
  mutate(unique_id = row_number())

# Define a color palette
color_palette <- c(
  "#1E3221", "#2A4731", "#395724", "#4B7038", "#5C884B",
  "#739D52", "#84A761", "#A2B9FC", "#8CB0F5", "#7583B7", 
  "#8A97CC", "#9FB5DF", "#DF4416", "#E2562C", "#E86D3D",
  "#EA8354", "#EB9C71", "#E1D83B", "#D66EB7", "#B05096", 
  "#C362A7", "#AB98C8", "#8A77AA", "#AFB4E2", "#C5D0F5",  
  "#4A2F6A", "#6A4E9C", "#BD3E6B", "#7C6A4B", "#9E7B5D"
)

# Create the facet wrapped plot
pro_HA_bar_plot <- ggplot(summary_df_pro, aes(x = Trimmed_Contrast, y = Number_of_Significant_proteins, fill = factor(unique_id))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(x = "Contrast", y = "Number of Significant Proteins") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 3, 3), "cm"),
    axis.title.y = element_text(size = 16, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)), 
    axis.title.x = element_text(size = 16, lineheight = .9, family = "sans", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1, lineheight = .9, family = "sans", colour = "black"),
    axis.text.y = element_text(size = 16, lineheight = .9, family = "sans", colour = "black"), 
    axis.title = element_text(size = 16, face = "bold", family = "sans", colour = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 16, lineheight = .9, family = "sans", colour = "black"),
    strip.text = element_text(size = 16, face = "bold", family = "sans", colour = "black")
  ) +
  facet_wrap(~ Category, scales = "free_x", ncol = 4)

print(pro_HA_bar_plot)

dir.create(OUTPUT_DIR, showWarnings = FALSE)
write.csv(summary_df_pro, file.path(OUTPUT_DIR, "figure_7a_pro_HA_bar_plot_data.csv"), row.names = FALSE)

path_output <- "figure_7a_pro_HA_bar_plot.pdf"
ggsave(path_output, pro_HA_bar_plot, width = 22, height = 8, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", path_output), pro_HA_bar_plot, width = 22, height = 8, dpi = 600, device = "tiff")
