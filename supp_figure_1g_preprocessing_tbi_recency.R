# =============================================================================
# Supplementary Figure 1g preprocessing – Biomarkers by time-since-mTBI category
# =============================================================================
#
# Description: Participants with an mTBI are split into tertiles of years
#   since their most recent mTBI (Recent / Intermediate / Remote), computed
#   once on the pooled mTBI (+) population so category boundaries are
#   identical whether looking at the pooled sample or a single group.
#   Plasma Abeta42/Abeta40 ratio and anti-Abeta1-42 IgG/IgA are compared
#   across these categories (Kruskal-Wallis omnibus, FDR-adjusted across
#   panels; Dunn post-hoc, Holm), both pooled and within each mTBI (+)
#   group. A continuous Spearman correlation between years since mTBI and
#   each biomarker is also computed as a robustness check.
#
#   Uses years since the most recent mTBI rather than the pipeline's
#   average_years_since_tbi: 53 of 72 mTBI (+) participants have 2+
#   recorded mTBIs, often decades apart, so averaging can blend injuries
#   into a value with no real injury near it. COVID-only participants (no
#   mTBI) are excluded, since time since mTBI is undefined for them.
#
# Prerequisites: supp_figure_1/supp_figure_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv;
#   supp_figure_1/supp_figure_1_input_data/tbi_most_recent.csv;
#   figure_1 output: figure_1i (N3PA Abeta ratio) in figure_1/figure_1_output_data/.
#
# Output: supp_figure_1_tbi_recency.pdf, supp_figure_1_tbi_recency_continuous.pdf;
#   data in supp_figure_1_output_data
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
input_dir  <- "supp_figure_1/supp_figure_1_input_data"
tif_dir    <- "supp_figure_1"
fig1_out   <- "figure_1/figure_1_output_data"

path_main <- file.path(input_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv")
path_recent_tbi <- file.path(input_dir, "tbi_most_recent.csv")
path_n3pa <- file.path(fig1_out, "figure_1i_N3PA_Abeta_plot_data.csv")
missing <- c(
  if (!file.exists(path_main)) "COAST_Study_Data_Clean_Age_Groups_add_dates.csv (in supp_figure_1/supp_figure_1_input_data)",
  if (!file.exists(path_recent_tbi)) "tbi_most_recent.csv (in supp_figure_1/supp_figure_1_input_data)",
  if (!file.exists(path_n3pa)) "figure_1i_N3PA_Abeta_plot_data.csv (run figure_1/figure_1i.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dunn.test)

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(path_main, stringsAsFactors = FALSE)
Tbi_Most_Recent <- read.csv(path_recent_tbi, stringsAsFactors = FALSE)
N3PA_SRX_Data_Clean_join <- read.csv(path_n3pa, stringsAsFactors = FALSE)

# mTBI (+) groups only; COVID-only (no mTBI) excluded (see header note).
tbi_groups <- c("COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)")

Corr_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::filter(qq_group %in% tbi_groups) %>%
  dplyr::select(participant_id, qq_group, vib_anti_abeta_1_42_igg_iga) %>%
  left_join(Tbi_Most_Recent %>% dplyr::select(participant_id, years_since_most_recent_tbi, n_tbi_events),
            by = "participant_id") %>%
  left_join(N3PA_SRX_Data_Clean_join %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R),
            by = c("participant_id", "qq_group")) %>%
  dplyr::filter(!is.na(years_since_most_recent_tbi))

cat("mTBI events per participant in this analysis:\n"); print(table(Corr_Data$n_tbi_events))

# Tertile cutpoints computed once on the pooled mTBI (+) population.
cuts <- quantile(Corr_Data$years_since_most_recent_tbi, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
Corr_Data$recency_category <- cut(Corr_Data$years_since_most_recent_tbi, breaks = cuts,
                                   include.lowest = TRUE, labels = c("Recent", "Intermediate", "Remote"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Corr_Data, file.path(output_dir, "supp_figure_1_tbi_recency_input_data.csv"), row.names = FALSE)

cat("Tertile cutpoints (years since mTBI):\n"); print(cuts)
cat("\nCategory counts:\n"); print(table(Corr_Data$recency_category, Corr_Data$qq_group))

# Long format: one row per participant x biomarker, plus a pooled "Overall" pseudo-group
biomarker_long <- Corr_Data %>%
  tidyr::pivot_longer(cols = c(Abeta42_Abeta40_Ratio_R, vib_anti_abeta_1_42_igg_iga),
                       names_to = "biomarker_name", values_to = "biomarker_value") %>%
  dplyr::mutate(biomarker_name = dplyr::recode(biomarker_name,
                                                Abeta42_Abeta40_Ratio_R = "Aβ42/Aβ40",
                                                vib_anti_abeta_1_42_igg_iga = "Anti-Aβ IgG/IgA")) %>%
  dplyr::filter(!is.na(biomarker_value), !is.na(recency_category))

plot_data <- dplyr::bind_rows(
  biomarker_long %>% dplyr::mutate(group_label = "Overall (mTBI +)"),
  biomarker_long %>% dplyr::mutate(group_label = qq_group)
)
plot_data$group_label <- factor(plot_data$group_label,
                                 levels = c("Overall (mTBI +)", tbi_groups))

# Kruskal-Wallis (omnibus) + Dunn post-hoc (Holm) per group_label x biomarker panel
compute_panel_stats <- function(df) {
  panels <- df %>% dplyr::distinct(group_label, biomarker_name)
  results <- lapply(seq_len(nrow(panels)), function(i) {
    g <- panels$group_label[i]; b <- panels$biomarker_name[i]
    sub <- df %>% dplyr::filter(group_label == g, biomarker_name == b)
    kw <- kruskal.test(biomarker_value ~ droplevels(recency_category), data = sub)
    dunn <- suppressWarnings(dunn.test(sub$biomarker_value, sub$recency_category, method = "holm", table = FALSE, list = FALSE))
    dunn_df <- data.frame(comparison = dunn$comparisons, Z = dunn$Z, p_value = dunn$P, p_holm = dunn$P.adjusted)
    dunn_df$group_label <- g
    dunn_df$biomarker_name <- b
    dunn_df$kruskal_p <- kw$p.value
    dunn_df$n <- nrow(sub)
    dunn_df
  })
  do.call(rbind, results)
}
panel_stats <- compute_panel_stats(plot_data)
# FDR across the omnibus Kruskal-Wallis tests (one per group_label x biomarker
# panel); Dunn post-hoc p_holm is already corrected within each panel.
kruskal_p_table <- unique(panel_stats[, c("group_label", "biomarker_name", "kruskal_p")])
kruskal_p_table$kruskal_p_fdr <- p.adjust(kruskal_p_table$kruskal_p, method = "fdr")
panel_stats <- merge(panel_stats, kruskal_p_table, by = c("group_label", "biomarker_name", "kruskal_p"))
panel_stats <- panel_stats[, c("group_label", "biomarker_name", "n", "kruskal_p", "kruskal_p_fdr", "comparison", "Z", "p_value", "p_holm")]
write.csv(panel_stats, file.path(output_dir, "supp_figure_1_tbi_recency_stats.csv"), row.names = FALSE)

Tbi_Recency_Plot <- ggplot(plot_data, aes(x = recency_category, y = biomarker_value)) +
  geom_boxplot(width = 0.5, fill = "white", color = "black", alpha = 0.7, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.15, height = 0, shape = 21, size = 2.5, stroke = 0.4, fill = "grey60") +
  facet_grid(biomarker_name ~ group_label, scales = "free_y") +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.5) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(size = 10, family = "Arial", face = "bold", colour = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10, family = "Arial", colour = "black"),
    axis.title.y = element_text(size = 13, family = "Arial", face = "bold", colour = "black"),
    plot.title = element_text(size = 16, family = "Arial", face = "bold")
  ) +
  labs(
    title = "Biomarkers by Time Since Most Recent mTBI (Tertile)",
    x = "Time Since Most Recent mTBI",
    y = "Biomarker Level"
  )

ggsave(file.path(tif_dir, "supp_figure_1_tbi_recency.pdf"), Tbi_Recency_Plot, width = 14, height = 8, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1_tbi_recency.pdf")), Tbi_Recency_Plot, width = 14, height = 8, dpi = 600, device = "tiff")

# -----------------------------------------------------------------------
# Continuous companion: Spearman correlation between years since mTBI and
# each biomarker (same panels as above), as a robustness check that
# categorizing time-since-exposure into tertiles did not change the
# conclusion.
# -----------------------------------------------------------------------
compute_cont_stats <- function(df) {
  panels <- df %>% dplyr::distinct(group_label, biomarker_name)
  results <- lapply(seq_len(nrow(panels)), function(i) {
    g <- panels$group_label[i]; b <- panels$biomarker_name[i]
    sub <- df %>% dplyr::filter(group_label == g, biomarker_name == b)
    ok <- stats::complete.cases(sub$years_since_most_recent_tbi, sub$biomarker_value)
    test <- suppressWarnings(stats::cor.test(sub$years_since_most_recent_tbi[ok], sub$biomarker_value[ok], method = "spearman"))
    data.frame(group_label = g, biomarker_name = b, n = sum(ok),
               rho = unname(test$estimate), p_value = test$p.value)
  })
  do.call(rbind, results)
}
cont_stats <- compute_cont_stats(plot_data)
cont_stats$p_fdr <- p.adjust(cont_stats$p_value, method = "fdr")
write.csv(cont_stats, file.path(output_dir, "supp_figure_1_tbi_recency_continuous_stats.csv"), row.names = FALSE)
cat("\nContinuous (Spearman) companion stats:\n"); print(cont_stats)

Tbi_Recency_Continuous_Plot <- ggplot(plot_data, aes(x = years_since_most_recent_tbi, y = biomarker_value)) +
  geom_point(shape = 21, size = 2.5, stroke = 0.4, fill = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  facet_grid(biomarker_name ~ group_label, scales = "free") +
  stat_cor(method = "spearman", size = 3.5) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, family = "Arial", face = "bold", colour = "black"),
    axis.text.x = element_text(size = 10, family = "Arial", colour = "black"),
    axis.text.y = element_text(size = 10, family = "Arial", colour = "black"),
    axis.title.y = element_text(size = 13, family = "Arial", face = "bold", colour = "black"),
    plot.title = element_text(size = 16, family = "Arial", face = "bold")
  ) +
  labs(
    title = "Biomarkers vs. Continuous Time Since Most Recent mTBI",
    x = "Years Since Most Recent mTBI",
    y = "Biomarker Level"
  )

ggsave(file.path(tif_dir, "supp_figure_1_tbi_recency_continuous.pdf"), Tbi_Recency_Continuous_Plot, width = 14, height = 8, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "supp_figure_1_tbi_recency_continuous.pdf")), Tbi_Recency_Continuous_Plot, width = 14, height = 8, dpi = 600, device = "tiff")
