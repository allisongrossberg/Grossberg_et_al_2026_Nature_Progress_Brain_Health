# =============================================================================
# Figure 1j – Anti-Aβ 1-42 IgG/IgA antibody levels by group
# =============================================================================
#
# Description: Box plot of anti-Aβ1-42 IgG/IgA by group.
#
#   Statistics: proportion above the control-group median, compared
#   pairwise across all six group pairs with Boschloo's exact test
#   (uncorrected for multiple comparisons; results for all six pairs are
#   saved to figure_1j_pairwise_stats.csv). Sidedness follows the same
#   pre-specified exposure-burden ordering as figure_1i.R (Control <
#   COVID-only/mTBI-only < Double+); note that the mTBI-only group's rate
#   on this biomarker (37.5% above the control median) happens to sit
#   below the control group's own rate (50%), so its one-sided test
#   against Control returns a large, non-significant p-value by
#   construction. Effect size (odds ratio, 95% CI) comes from a two-sided
#   Fisher's exact test on the same 2x2 table. The caption highlights the
#   significant pair, COVID-19 (+) mTBI (+) vs. COVID-19 (-) mTBI (+)
#   (p = 0.033).
#
# Prerequisites: figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: figure_1j_Anti_Abeta_IgG_IgA.pdf; data in figure_1_output_data
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
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(
  file.path(output_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"),
  stringsAsFactors = FALSE
)

library(dplyr)
library(plotrix)
library(exact2x2)

#######

Vib_Figures_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% dplyr::select(qq_group, participant_id, vib_anti_abeta_1_42_igg_iga, vib_anti_abeta_1_42_igm)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Vib_Figures_Data %>% dplyr::select(participant_id, qq_group, vib_anti_abeta_1_42_igg_iga, vib_anti_abeta_1_42_igm),
          file.path(output_dir, "figure_1j_Vib_Figures_Data.csv"), row.names = FALSE)

library(tidyverse)
library(ggsignif)
library(cowplot)

Vib_Figures_Data_long <- Vib_Figures_Data %>%
  pivot_longer(
    cols = c(vib_anti_abeta_1_42_igg_iga, vib_anti_abeta_1_42_igm),
    names_to = "label",
    values_to = "value"
  )

Vib_Figures_Data_long <- Vib_Figures_Data_long %>%
  mutate(label = case_when(
    label == "vib_anti_abeta_1_42_igg_iga" ~ "IgG/IgA",
    label == "vib_anti_abeta_1_42_igm" ~ "IgM",
    TRUE ~ label))

aggregate(vib_anti_abeta_1_42_igg_iga ~ qq_group,
          data = COAST_Study_Data_Clean_Age_Groups_add_dates,
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = std.error(x)))

# Headline statistic: pairwise Boschloo's exact test (see header).
median_value_control <- Vib_Figures_Data %>%
  filter(qq_group == "COVID-19 (-) mTBI (-)") %>%
  summarize(median_value = median(vib_anti_abeta_1_42_igg_iga, na.rm = TRUE)) %>%
  pull(median_value)
cat("\nControl-group median (threshold):", median_value_control, "\n")

Vib_Figures_Data$above_median <- Vib_Figures_Data$vib_anti_abeta_1_42_igg_iga > median_value_control

group_counts_above <- Vib_Figures_Data %>%
  dplyr::filter(!is.na(above_median)) %>%
  dplyr::group_by(qq_group) %>%
  dplyr::summarise(n = dplyr::n(), x_above = sum(above_median), .groups = "drop")
print(group_counts_above)

format_pvalue <- function(p) {
  if (p < 0.0001) return("p < 0.0001")
  return(paste("p =", formatC(p, format = "f", digits = 4)))
}
format_ci <- function(lower, upper) {
  upper_str <- if (is.infinite(upper)) "∞" else sprintf("%.2f", upper)
  sprintf("[%.2f, %s]", lower, upper_str)
}

# Exposure-burden rank: Control (0) < COVID-only/mTBI-only (1, tied) <
# Double+ (2). Used below to assign one-sided vs. two-sided per pair.
group_rank <- c(
  "COVID-19 (-) mTBI (-)" = 0,
  "COVID-19 (+) mTBI (-)" = 1,
  "COVID-19 (-) mTBI (+)" = 1,
  "COVID-19 (+) mTBI (+)" = 2
)

compute_pairwise_boschloo <- function(group_counts, count_col, group_rank) {
  pairs <- utils::combn(group_counts$qq_group, 2, simplify = FALSE)
  results <- lapply(pairs, function(p) {
    rank1 <- group_rank[[p[1]]]; rank2 <- group_rank[[p[2]]]
    if (rank1 == rank2) {
      alt <- "two.sided"
    } else {
      if (rank2 > rank1) p <- rev(p)  # put the higher-rank group first
      alt <- "less"
    }
    r1 <- group_counts[group_counts$qq_group == p[1], ]
    r2 <- group_counts[group_counts$qq_group == p[2], ]
    x1 <- r1[[count_col]]; n1 <- r1$n; x2 <- r2[[count_col]]; n2 <- r2$n
    # boschloo()'s "alternative" refers to the ratio p2(1-p1)/[p1(1-p2)];
    # "less" (ratio < 1) corresponds to p1 > p2, i.e. group_1 (the
    # higher-rank group, ordered first above) has the higher rate.
    bosch <- boschloo(x1, n1, x2, n2, alternative = alt)
    fish  <- fisher.test(matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE))
    data.frame(group_1 = p[1], group_2 = p[2], x1 = x1, n1 = n1, x2 = x2, n2 = n2,
               alternative = alt,
               boschloo_p = bosch$p.value, fisher_p_twosided = fish$p.value,
               or_twosided = unname(fish$estimate),
               or_ci_lower_twosided = fish$conf.int[1], or_ci_upper_twosided = fish$conf.int[2])
  })
  out <- do.call(rbind, results)
  out[order(out$boschloo_p), ]
}

pairwise_stats <- compute_pairwise_boschloo(group_counts_above, "x_above", group_rank)
cat("\nAll 6 pairwise comparisons (proportion above control median), sorted by p-value:\n")
print(pairwise_stats, digits = 4, row.names = FALSE)

write.csv(pairwise_stats, file.path(output_dir, "figure_1j_pairwise_stats.csv"), row.names = FALSE)

# Highlighted comparison for the caption: Double+ vs. mTBI-only.
highlight <- pairwise_stats[
  (pairwise_stats$group_1 == "COVID-19 (+) mTBI (+)" & pairwise_stats$group_2 == "COVID-19 (-) mTBI (+)") |
  (pairwise_stats$group_2 == "COVID-19 (+) mTBI (+)" & pairwise_stats$group_1 == "COVID-19 (-) mTBI (+)"), ]
# Orient so Double+ (the higher-rate group) is reported first.
if (highlight$group_1 == "COVID-19 (+) mTBI (+)") {
  x_focal <- highlight$x1; n_focal <- highlight$n1; x_other <- highlight$x2; n_other <- highlight$n2
} else {
  x_focal <- highlight$x2; n_focal <- highlight$n2; x_other <- highlight$x1; n_other <- highlight$n1
}
or_estimate <- highlight$or_twosided
or_ci_lower <- highlight$or_ci_lower_twosided
or_ci_upper <- highlight$or_ci_upper_twosided
boschloo_p  <- highlight$boschloo_p

cat(sprintf("\nHighlighted pair - Double+ vs mTBI-only: %d/%d vs %d/%d above median, one-sided Boschloo p = %.4f\n",
            x_focal, n_focal, x_other, n_other, boschloo_p))

caption_text <- sprintf(
  "Boschloo's exact test (one-sided): %d/%d vs %d/%d above median; OR = %.2f, 95%% CI %s; %s",
  x_focal, n_focal, x_other, n_other, or_estimate, format_ci(or_ci_lower, or_ci_upper),
  format_pvalue(boschloo_p)
)

stats_out <- data.frame(
  highlighted_pair = "COVID-19 (+) mTBI (+) vs COVID-19 (-) mTBI (+)",
  threshold_type = "control-group median", threshold_value = median_value_control,
  x_focal = x_focal, n_focal = n_focal, x_other = x_other, n_other = n_other,
  boschloo_p_twosided = boschloo_p,
  or_twosided = or_estimate, or_ci_lower_twosided = or_ci_lower, or_ci_upper_twosided = or_ci_upper,
  note = "See figure_1j_pairwise_stats.csv for all 6 uncorrected pairwise comparisons."
)
write.csv(stats_out, file.path(output_dir, "figure_1j_stats.csv"), row.names = FALSE)

Vib_Figures_Data_long$qq_group <- factor(Vib_Figures_Data_long$qq_group, levels = c("COVID-19 (-) mTBI (-)", "COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))
plot_igg_iga <- ggplot(Vib_Figures_Data_long %>% filter(label == "IgG/IgA"), aes(x = qq_group, y = value)) +
  # Control-group median: the threshold the Boschloo's exact test above is computed against.
  geom_hline(yintercept = median_value_control, linetype = "dotted", color = "#DD4726", linewidth = 3) +
  geom_boxplot(width = 0.5, position = position_dodge(0.8), fill = "white", color = "black", alpha = 0.7, outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(aes(fill = qq_group), width = 0.2, height = 0, shape = 21, size = 4, stroke = 0.5) +
  scale_fill_brewer(palette = "Greys") +
  theme_minimal() +
  theme(legend.title=element_blank(),legend.position="none",
      plot.margin = unit(c(1,1,1,1), "cm"),
      panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black"),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.text.x = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
      axis.text.y = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
      legend.text = element_text(size = 18, lineheight = .9, family = "Arial", face = "bold", colour = "black"),
      plot.caption = element_text(size = 14, face = "bold", margin = margin(t = 20), family = "Arial"))+
  labs(
    y = "Anti-Aβ1-42 (IgG/IgA)",
    x = "Group",
    caption = caption_text
  ) +
  scale_x_discrete(labels = c("COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)",
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)"))+
  coord_cartesian(clip = "off") +
  # Explicit breaks so this axis reads identically to supp_figure_1d.R (same variable, same pooled data).
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25))

ggsave(file.path(tif_dir, "figure_1j_Anti_Abeta_IgG_IgA.pdf"), plot_igg_iga, width = 12, height = 10, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_1j_Anti_Abeta_IgG_IgA.pdf")), plot_igg_iga, width = 12, height = 10, dpi = 600, device = "tiff")
