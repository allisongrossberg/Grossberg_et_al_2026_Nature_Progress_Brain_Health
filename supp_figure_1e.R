# =============================================================================
# Supplementary Figure 1e – Correlation heatmap: symptoms, clinical outcomes,
#   and biomarkers (symptomatic groups combined)
# =============================================================================
#
# Description: Spearman rank correlations between total neuro/psych symptom
#   severity/frequency scores and PHQ-8, EQ-5D, WAI-2, Neuro-QOL, plasma
#   Aβ42/Aβ40 ratio, and anti-Aβ1-42 IgG/IgA, pooled across the three
#   symptomatic groups (COVID-19+mTBI-, COVID-19-mTBI+, COVID-19+mTBI+). The
#   double-negative control group is excluded (zero variance in symptom
#   severity/frequency, so correlation is undefined for that group).
#
# Prerequisites: supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv;
#   figure_1 outputs: figure_1f (NeuroQOL), figure_1i (N3PA), figure_1h (PHQ8),
#   figure_1e (WAI) in figure_1/figure_1_output_data/.
#
# Output: supp_figure_1e_Correlation_Heatmap.pdf; supp_figure_1_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
if (basename(getwd()) == "supp_figure_1") setwd("..")
output_dir <- "supp_figure_1/supp_figure_1_output_data"
input_dir  <- "supp_figure_1/supp_figure_1_input_data"
tif_dir    <- "supp_figure_1"
fig1_out   <- "figure_1/figure_1_output_data"

path_main <- file.path(input_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv")
path_nq   <- file.path(fig1_out, "figure_1f_NeuroQOL_plot_data.csv")
path_n3pa <- file.path(fig1_out, "figure_1i_N3PA_Abeta_plot_data.csv")
path_phq8 <- file.path(fig1_out, "figure_1h_PHQ8_plot_data.csv")
path_wai  <- file.path(fig1_out, "figure_1e_WAI_plot_data.csv")
missing <- c(
  if (!file.exists(path_main)) "COAST_Study_Data_Clean_Age_Groups_add_dates.csv (in supp_figure_1_input_data)",
  if (!file.exists(path_nq))   "figure_1f_NeuroQOL_plot_data.csv (run figure_1/figure_1f.R)",
  if (!file.exists(path_n3pa)) "figure_1i_N3PA_Abeta_plot_data.csv (run figure_1/figure_1i.R)",
  if (!file.exists(path_phq8)) "figure_1h_PHQ8_plot_data.csv (run figure_1/figure_1h.R)",
  if (!file.exists(path_wai))  "figure_1e_WAI_plot_data.csv (run figure_1/figure_1e.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(path_main, stringsAsFactors = FALSE)
NeuroQOL_Data_join       <- read.csv(path_nq, stringsAsFactors = FALSE)
N3PA_SRX_Data_Clean_join <- read.csv(path_n3pa, stringsAsFactors = FALSE)
PHQ8_plot_data           <- read.csv(path_phq8, stringsAsFactors = FALSE)
WAI_plot_data            <- read.csv(path_wai, stringsAsFactors = FALSE)

# Symptomatic groups only; double-negative control excluded (see header note).
group_levels <- c("COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)")

Corr_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::filter(qq_group %in% group_levels) %>%
  dplyr::select(participant_id, qq_group,
                recode_overall_neuro_psych_sev_score, recode_overall_neuro_psych_freq_score,
                qq_eq5d_index_score, vib_anti_abeta_1_42_igg_iga) %>%
  left_join(NeuroQOL_Data_join %>% dplyr::select(participant_id, qq_group, Neuro_QOL_TScore),
            by = c("participant_id", "qq_group")) %>%
  left_join(N3PA_SRX_Data_Clean_join %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R),
            by = c("participant_id", "qq_group")) %>%
  left_join(PHQ8_plot_data %>% dplyr::select(participant_id, qq_group, qq_phq8_average_score),
            by = c("participant_id", "qq_group")) %>%
  left_join(WAI_plot_data %>% dplyr::select(participant_id, qq_group, qq_wai_2),
            by = c("participant_id", "qq_group"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Corr_Data, file.path(output_dir, "supp_figure_1e_input_data.csv"), row.names = FALSE)

# Variables of interest and display labels
vars <- c("recode_overall_neuro_psych_sev_score", "recode_overall_neuro_psych_freq_score",
          "qq_phq8_average_score", "qq_eq5d_index_score", "qq_wai_2",
          "Neuro_QOL_TScore", "Abeta42_Abeta40_Ratio_R", "vib_anti_abeta_1_42_igg_iga")
labels <- c("Total Neuro/Psych\nSeverity", "Total Neuro/Psych\nFrequency",
            "PHQ-8", "EQ-5D", "WAI-2", "Neuro-QOL",
            "Aβ42/Aβ40", "Anti-Aβ\nIgG/IgA")

# Pairwise Spearman correlation with p-values, robust to small/unequal n per pair.
# (cor.test requires n >= 3; pairs with fewer complete observations are left NA.)
compute_cor_stats <- function(df, vars) {
  n_vars <- length(vars)
  cor_mat <- matrix(NA_real_, n_vars, n_vars, dimnames = list(vars, vars))
  p_mat   <- matrix(NA_real_, n_vars, n_vars, dimnames = list(vars, vars))
  n_mat   <- matrix(NA_real_, n_vars, n_vars, dimnames = list(vars, vars))
  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {
      x <- df[[vars[i]]]; y <- df[[vars[j]]]
      ok <- stats::complete.cases(x, y)
      n_pairs <- sum(ok)
      n_mat[i, j] <- n_pairs
      if (i == j) {
        cor_mat[i, j] <- 1
        p_mat[i, j] <- 0
      } else if (n_pairs >= 3) {
        test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
        cor_mat[i, j] <- unname(test$estimate)
        p_mat[i, j] <- test$p.value
      }
    }
  }
  list(cor = cor_mat, p = p_mat, n = n_mat)
}

# Long-format pairwise table (upper triangle only) with FDR-adjusted p-values
build_pairwise_table <- function(stats_list, vars) {
  idx <- utils::combn(vars, 2)
  out <- lapply(seq_len(ncol(idx)), function(k) {
    v1 <- idx[1, k]; v2 <- idx[2, k]
    data.frame(variable_1 = v1, variable_2 = v2,
               rho = stats_list$cor[v1, v2], p_value = stats_list$p[v1, v2],
               n = stats_list$n[v1, v2], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  out$p_fdr <- p.adjust(out$p_value, method = "fdr")
  out[order(out$p_value), ]
}

overall_stats <- compute_cor_stats(Corr_Data, vars)
overall_pairwise <- build_pairwise_table(overall_stats, vars)
write.csv(overall_pairwise, file.path(output_dir, "supp_figure_1e_pairwise_stats.csv"), row.names = FALSE)

n_participants <- length(unique(Corr_Data$participant_id))
n_participants

# Display matrices with pretty labels for plotting
cor_mat_disp <- overall_stats$cor
p_mat_disp   <- overall_stats$p
dimnames(cor_mat_disp) <- list(labels, labels)
dimnames(p_mat_disp)   <- list(labels, labels)

# Coefficient labels, with a p-value-tiered "*" superscript directly on the
# number (via plotmath) -- so the flag sits right on the coefficient itself,
# not off in the corner of the tile. Tiers: * p<0.05, ** p<0.01, *** p<0.001.
add_coefficient_labels <- function(plot, cor_mat, p_mat, digits = 2, size = 4) {
  lower_cor <- cor_mat
  lower_cor[upper.tri(lower_cor, diag = TRUE)] <- NA
  lower_p <- p_mat
  lower_p[upper.tri(lower_p, diag = TRUE)] <- NA

  cor_long <- reshape2::melt(lower_cor, na.rm = TRUE)
  p_long   <- reshape2::melt(lower_p, na.rm = TRUE)
  lab_df <- cor_long
  lab_df$p <- p_long$value
  formatted <- as.character(round(lab_df$value, digits))
  stars <- ifelse(lab_df$p < 0.001, "***",
            ifelse(lab_df$p < 0.01,  "**",
             ifelse(lab_df$p < 0.05,  "*", "")))
  lab_df$label <- ifelse(
    stars == "",
    paste0('"', formatted, '"'),
    paste0('"', formatted, '"^"', stars, '"')
  )

  plot + ggplot2::geom_text(
    data = lab_df,
    mapping = ggplot2::aes(x = Var1, y = Var2, label = label),
    parse = TRUE, size = size
  )
}

Correlelogram_Overall_Plot <- ggcorrplot(
  cor_mat_disp,
  type = "lower",
  outline.color = "white",
  ggtheme = ggplot2::theme_minimal,
  colors = c("#7583B7", "white", "#DD4726"),
  lab = FALSE,
  legend.title = "Spearman's ρ"
) +
  theme(
    axis.text.x = element_text(size = 13, family = "Arial", face = "bold", colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13, family = "Arial", face = "bold", colour = "black"),
    legend.title = element_text(size = 12, family = "Arial", face = "bold"),
    legend.text = element_text(size = 11, family = "Arial")
  )
Correlelogram_Overall_Plot <- add_coefficient_labels(Correlelogram_Overall_Plot, cor_mat_disp, p_mat_disp, size = 6)

ggsave(file.path(tif_dir, "supp_figure_1e_Correlation_Heatmap.pdf"), Correlelogram_Overall_Plot,
       width = 10, height = 9, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1e_Correlation_Heatmap.pdf")), Correlelogram_Overall_Plot,
       width = 10, height = 9, dpi = 600, device = "tiff")
