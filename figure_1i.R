# =============================================================================
# Figure 1i – Plasma Aβ42/Aβ40 ratio by group
# =============================================================================
#
# Description: Box plot of plasma Aβ42/Aβ40 ratio (N3PA) across four groups.
#   Reference line at 0.04 (literature-derived amyloid-abnormality cutoff).
#   Run preprocessing first.
#
#   Statistics: proportion below 0.04 (amyloid-abnormality threshold),
#   compared pairwise across all six group pairs with Boschloo's exact test
#   (uncorrected for multiple comparisons; results for all six pairs are
#   saved to figure_1i_pairwise_stats.csv). Sidedness follows a
#   pre-specified exposure-burden ordering, Control < COVID-only/mTBI-only
#   < Double+: pairs with different rank are tested one-sided (higher rank
#   > lower rank); the one tied-rank pair (COVID-only vs. mTBI-only) is
#   two-sided. Effect size (odds ratio, 95% CI) comes from a two-sided
#   Fisher's exact test on the same 2x2 table. The caption highlights the
#   significant pair, COVID-19 (+) mTBI (+) vs. control (one-sided
#   p = 0.0213).
#
# Prerequisites: figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv,
#   figure_1_input_data/N3PA_SRX_Data_7_30_24.xlsx (run preprocessing step 1 & 2 first).
#
# Output: figure_1i_Abeta_Ratio.pdf; data in figure_1_output_data
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
output_dir <- if (dir.exists("figure_1_revised")) "figure_1_revised/figure_1_output_data" else "figure_1_output_data"
input_dir  <- if (dir.exists("figure_1_revised")) "figure_1_revised/figure_1_input_data" else "figure_1_input_data"
tif_dir    <- if (dir.exists("figure_1_revised")) "figure_1_revised" else "."
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(
  file.path(output_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"),
  stringsAsFactors = FALSE
)
join_and_select <- function(df_A, df_B, columns_to_select) {
  if (!("participant_id" %in% names(df_A)) || !("participant_id" %in% names(df_B))) stop("Both dataframes must have 'participant_id'")
  if (!all(columns_to_select %in% names(df_A))) stop("columns_to_select must be in df_A")
  df_A %>% dplyr::select(participant_id, dplyr::all_of(columns_to_select)) %>% dplyr::right_join(df_B, by = "participant_id")
}

library(ggplot2)
library(dplyr)
library(plotrix)
library(readxl)
library(exact2x2)

# Use N3PA CSV
n3pa_csv  <- file.path(input_dir, "N3PA_SRX_Data_7_30_24.csv")
n3pa_xlsx <- file.path(input_dir, "N3PA_SRX_Data_7_30_24.xlsx")
if (file.exists(n3pa_csv)) {
  N3PA_SRX_Data <- read.csv(n3pa_csv, stringsAsFactors = FALSE)
} else {
  N3PA_SRX_Data <- read_excel(n3pa_xlsx)
}
N3PA_SRX_Data_Clean <- N3PA_SRX_Data %>% filter((N3PA_Status == "Complete"))

N3PA_SRX_Data_Clean_join <- join_and_select(COAST_Study_Data_Clean_Age_Groups_add_dates,
                                            N3PA_SRX_Data_Clean,
                                            c("qq_biological_sex",
                                              "age_years",
                                              "qq_group",
                                              "recode_overall_neuro_psych_sev_score"))


N3PA_SRX_Data_Clean_join$Abeta42_Abeta40_Ratio_R <- N3PA_SRX_Data_Clean_join$N3PA_Abeta42_Concentration / N3PA_SRX_Data_Clean_join$N3PA_Abeta40_Concentration

N3PA_SRX_Data_Clean_join <- N3PA_SRX_Data_Clean_join %>% filter(!is.na(qq_group))

N3PA_SRX_Data_Clean_join <- N3PA_SRX_Data_Clean_join %>% filter(!is.na(Abeta42_Abeta40_Ratio_R))

N3PA_plot_data <- N3PA_SRX_Data_Clean_join %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(N3PA_plot_data %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R), file.path(output_dir, "figure_1i_N3PA_Abeta_plot_data.csv"), row.names = FALSE)

summary_data <- N3PA_SRX_Data_Clean_join %>%
  group_by(qq_group) %>%
  dplyr::summarise(mean_y = mean(Abeta42_Abeta40_Ratio_R),
                   sd_y = sd(Abeta42_Abeta40_Ratio_R),
                   se_y = sd(Abeta42_Abeta40_Ratio_R) / sqrt(n()))

aggregate(Abeta42_Abeta40_Ratio_R ~ qq_group,
          data = N3PA_SRX_Data_Clean_join,
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = std.error(x)))

N3PA_SRX_Data_Clean_join$Abeta_Group <- ifelse(N3PA_SRX_Data_Clean_join$Abeta42_Abeta40_Ratio_R < 0.04, "Below_0.04", "Above_0.04")

# Headline statistic: pairwise Boschloo's exact test, proportion below
# 0.04 (see header for sidedness rule).
group_counts_below <- N3PA_SRX_Data_Clean_join %>%
  dplyr::group_by(qq_group) %>%
  dplyr::summarise(n = dplyr::n(), x_below = sum(Abeta42_Abeta40_Ratio_R < 0.04), .groups = "drop")
print(group_counts_below)

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
    # boschloo()'s alternative refers to the ratio p2(1-p1)/[p1(1-p2)];
    # "less" (ratio < 1) means p1 > p2, i.e. group_1 has the higher rate.
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

pairwise_stats <- compute_pairwise_boschloo(group_counts_below, "x_below", group_rank)
cat("\nAll 6 pairwise comparisons (proportion below 0.04), sorted by p-value:\n")
print(pairwise_stats, digits = 4, row.names = FALSE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(pairwise_stats, file.path(output_dir, "figure_1i_pairwise_stats.csv"), row.names = FALSE)

# Caption highlights Double+ vs. Control (significant) and COVID-only vs.
# Control (p = 0.0522).
extract_pair <- function(g1, g2) {
  pairwise_stats[pairwise_stats$group_1 == g1 & pairwise_stats$group_2 == g2, ]
}
double_vs_control <- extract_pair("COVID-19 (+) mTBI (+)", "COVID-19 (-) mTBI (-)")
covid_vs_control   <- extract_pair("COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (-)")

cat(sprintf("\nHighlighted pair - Double+ vs Control: %d/%d vs %d/%d below 0.04, one-sided Boschloo p = %.4f\n",
            double_vs_control$x1, double_vs_control$n1, double_vs_control$x2, double_vs_control$n2,
            double_vs_control$boschloo_p))
cat(sprintf("Highlighted pair - COVID-only vs Control: %d/%d vs %d/%d below 0.04, one-sided Boschloo p = %.4f\n",
            covid_vs_control$x1, covid_vs_control$n1, covid_vs_control$x2, covid_vs_control$n2,
            covid_vs_control$boschloo_p))

# strwrap() balances the two caption lines to near-equal length for clean
# right alignment.
caption_full <- sprintf(
  "Boschloo's exact test (one-sided), vs. control (%d/%d): COVID-19 (+) mTBI (+) %d/%d, OR = %.2f, 95%% CI %s, %s; COVID-19 (+) mTBI (-) %d/%d, OR = %.2f, 95%% CI %s, %s",
  double_vs_control$x2, double_vs_control$n2,
  double_vs_control$x1, double_vs_control$n1,
  double_vs_control$or_twosided, format_ci(double_vs_control$or_ci_lower_twosided, double_vs_control$or_ci_upper_twosided),
  format_pvalue(double_vs_control$boschloo_p),
  covid_vs_control$x1, covid_vs_control$n1,
  covid_vs_control$or_twosided, format_ci(covid_vs_control$or_ci_lower_twosided, covid_vs_control$or_ci_upper_twosided),
  format_pvalue(covid_vs_control$boschloo_p)
)
caption_text <- paste(strwrap(caption_full, width = 105), collapse = "\n")

# Significant pair (Double+ vs. Control), for the stats CSV below.
x_focal <- double_vs_control$x1; n_focal <- double_vs_control$n1
x_other <- double_vs_control$x2; n_other <- double_vs_control$n2
or_estimate <- double_vs_control$or_twosided
or_ci_lower <- double_vs_control$or_ci_lower_twosided
or_ci_upper <- double_vs_control$or_ci_upper_twosided
boschloo_p  <- double_vs_control$boschloo_p

library(ggbreak)
library(ggh4x)
library(gg.gap)

y1_end <- 0
y2_start <- 0.02

N3PA_SRX_Data_Clean_join$qq_group <- factor(N3PA_SRX_Data_Clean_join$qq_group, levels = c("COVID-19 (-) mTBI (-)", "COVID-19 (+) mTBI (-)", "COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))

# Upper bound set from the data (not hardcoded) so all points are shown;
# a fixed 0-0.06 range would clip 18 of 97 participants whose ratio
# exceeds 0.06, distorting the apparent below/above-0.04 balance.
max_ratio_value <- max(N3PA_SRX_Data_Clean_join$Abeta42_Abeta40_Ratio_R, na.rm = TRUE)
y_upper <- ceiling(max_ratio_value / 0.01) * 0.01

Abeta_Ratio_Final_Plot <- ggplot(N3PA_SRX_Data_Clean_join, aes(x = qq_group, y = Abeta42_Abeta40_Ratio_R, fill = factor(qq_group)))+
  geom_hline(yintercept = 0.04, linetype = "dotted", color = "#DD4726", linewidth = 3)+
  geom_boxplot(position = position_dodge(0.8), width = 0.5, alpha = 0.7, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.8)+
  geom_jitter(aes(fill = qq_group), width = 0.2, shape = 21, size=4, stroke = 0.5)+
  scale_fill_brewer(palette = "Greys") +
  theme_minimal()+
  scale_y_continuous(limits = c(0, y_upper), breaks = seq(0, y_upper, by = 0.01))+
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
  labs(y = "Plasma Aβ42/Aβ40 Ratio",
       x = "Group",
       caption = caption_text)+
  scale_x_discrete(labels = c("COVID-19 (-)\nmTBI (-)", "COVID-19 (+)\nmTBI (-)",
                              "COVID-19 (-)\nmTBI (+)", "COVID-19 (+)\nmTBI (+)")) +
  coord_cartesian(clip = "off")


participant_count <- length(unique(N3PA_SRX_Data_Clean_join$participant_id))
participant_count

# Group-wise percent below threshold, for reference.
group_counts <- N3PA_SRX_Data_Clean_join %>%
  group_by(qq_group) %>%
  dplyr::summarise(total_individuals = n())
below_threshold_counts <- N3PA_SRX_Data_Clean_join %>%
  filter(Abeta42_Abeta40_Ratio_R < 0.04) %>%
  group_by(qq_group) %>%
  dplyr::summarise(below_threshold_individuals = n())
percentage_df <- group_counts %>%
  left_join(below_threshold_counts, by = "qq_group") %>%
  mutate(Percentage = (below_threshold_individuals / total_individuals) * 100) %>%
  select(Group = qq_group, Percentage)
print(percentage_df)

stats_out <- data.frame(
  highlighted_pair = "COVID-19 (+) mTBI (+) vs COVID-19 (-) mTBI (-)",
  threshold = 0.04, x_focal = x_focal, n_focal = n_focal, x_other = x_other, n_other = n_other,
  boschloo_p_twosided = boschloo_p,
  or_twosided = or_estimate, or_ci_lower_twosided = or_ci_lower, or_ci_upper_twosided = or_ci_upper,
  note = "See figure_1i_pairwise_stats.csv for all 6 uncorrected pairwise comparisons."
)
write.csv(stats_out, file.path(output_dir, "figure_1i_stats.csv"), row.names = FALSE)

ggsave(file.path(tif_dir, "figure_1i_Abeta_Ratio.pdf"), Abeta_Ratio_Final_Plot, width = 12, height = 10, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_1i_Abeta_Ratio.pdf")), Abeta_Ratio_Final_Plot, width = 12, height = 10, dpi = 600, device = "tiff")
