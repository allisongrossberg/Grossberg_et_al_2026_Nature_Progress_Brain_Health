# =============================================================================
# Supplementary Table 4 – Clinical and biomarker outcomes by COVID-19, mTBI, and sex
# =============================================================================
#
# Description: Table of outcome and biomarker variables (mean [SEM]) by study
#   group and sex. Wilcoxon rank-sum (Female vs Male within each group).
#   Requires joined data from main COAST file and figure 1f/1i/1j outputs.
#
# Prerequisites: supp_table_4_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv;
#   figure_1 outputs: figure_1f (NeuroQOL), figure_1i (N3PA), figure_1j (Vib),
#   figure_1h (PHQ8), figure_1e (WAI) in figure_1/figure_1_output_data/.
#
# Output: supp_table_4/supp_table_4.docx; supp_table_4_output_data/supp_table_4_export.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) %in% paste0("supp_table_", 1:5)) setwd("..")
in_dir <- "supp_table_4/supp_table_4_input_data"
fig1_out <- "figure_1/figure_1_output_data"
path_main <- file.path(in_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv")
path_nq   <- file.path(fig1_out, "figure_1f_NeuroQOL_plot_data.csv")
path_n3pa <- file.path(fig1_out, "figure_1i_N3PA_Abeta_plot_data.csv")
path_vib  <- file.path(fig1_out, "figure_1j_Vib_Figures_Data.csv")
path_phq8 <- file.path(fig1_out, "figure_1h_PHQ8_plot_data.csv")
path_wai  <- file.path(fig1_out, "figure_1e_WAI_plot_data.csv")
missing <- c(
  if (!file.exists(path_main)) "COAST_Study_Data_Clean_Age_Groups_add_dates.csv (in supp_table_4_input_data)",
  if (!file.exists(path_nq))   "figure_1f_NeuroQOL_plot_data.csv (run figure_1/figure_1f.R)",
  if (!file.exists(path_n3pa)) "figure_1i_N3PA_Abeta_plot_data.csv (run figure_1/figure_1i.R)",
  if (!file.exists(path_vib))  "figure_1j_Vib_Figures_Data.csv (run figure_1/figure_1j.R)",
  if (!file.exists(path_phq8)) "figure_1h_PHQ8_plot_data.csv (run figure_1/figure_1h.R)",
  if (!file.exists(path_wai))  "figure_1e_WAI_plot_data.csv (run figure_1/figure_1e.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(path_main, stringsAsFactors = FALSE)
NeuroQOL_Data_join <- read.csv(path_nq, stringsAsFactors = FALSE)
N3PA_SRX_Data_Clean_join <- read.csv(path_n3pa, stringsAsFactors = FALSE)
Vib_Figures_Data <- read.csv(path_vib, stringsAsFactors = FALSE)
PHQ8_plot_data <- read.csv(path_phq8, stringsAsFactors = FALSE)
WAI_plot_data <- read.csv(path_wai, stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(flextable)
library(officer)

# Join outcome and biomarker data (COAST + figure_1 outputs)
Table_4_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  dplyr::select(participant_id, qq_group, qq_biological_sex,
                recode_overall_neuro_psych_sev_score, recode_overall_neuro_psych_freq_score,
                qq_eq5d_index_score) %>%
  left_join(NeuroQOL_Data_join %>% dplyr::select(participant_id, qq_group, Neuro_QOL_TScore),
            by = c("participant_id", "qq_group")) %>%
  left_join(N3PA_SRX_Data_Clean_join %>% dplyr::select(participant_id, qq_group, Abeta42_Abeta40_Ratio_R),
            by = c("participant_id", "qq_group")) %>%
  left_join(Vib_Figures_Data %>% dplyr::select(participant_id, qq_group, vib_anti_abeta_1_42_igg_iga),
            by = c("participant_id", "qq_group")) %>%
  left_join(PHQ8_plot_data %>% dplyr::select(participant_id, qq_group, qq_phq8_average_score),
            by = c("participant_id", "qq_group")) %>%
  left_join(WAI_plot_data %>% dplyr::select(participant_id, qq_group, qq_wai_2),
            by = c("participant_id", "qq_group"))

# Outcome variables in table order (name = display label)
outcome_vars <- c(
  recode_overall_neuro_psych_sev_score = "Total Neurological/Psychological Symptom Severity Score",
  recode_overall_neuro_psych_freq_score = "Total Neurological/Psychological Symptom Frequency Score",
  qq_phq8_average_score = "PHQ-8 Score",
  qq_eq5d_index_score = "EQ-5D Score",
  Neuro_QOL_TScore = "Neuro-QOL T-Score",
  qq_wai_2 = "WAI Score",
  Abeta42_Abeta40_Ratio_R = "Aβ42/Aβ40",
  vib_anti_abeta_1_42_igg_iga = "Anti-Aβ IgG/IgA"
)

group_levels <- c(
  "COVID-19 (-) mTBI (-)",
  "COVID-19 (+) mTBI (-)",
  "COVID-19 (-) mTBI (+)",
  "COVID-19 (+) mTBI (+)"
)

Table_4_Data$qq_group <- factor(Table_4_Data$qq_group, levels = group_levels)

# Compute mean (SEM) and Wilcoxon p-value for Female vs Male within each qq_group
mean_sem <- function(x) {
  x <- na.omit(as.numeric(x))
  n <- length(x)
  if (n == 0) return(list(mean = NA, sem = NA, n = 0))
  m <- mean(x)
  s <- if (n > 1) sd(x) / sqrt(n) else 0
  list(mean = m, sem = s, n = n)
}

format_mean_sem <- function(m, sem, digits = 2) {
  if (is.na(m) || is.na(sem)) return("N/A")
  paste0(sprintf(paste0("%.", digits, "f"), m), " (", sprintf(paste0("%.", digits, "f"), sem), ")")
}

build_block <- function(data, var) {
  data <- data %>% filter(!is.na(!!sym(var)))
  female <- data %>% filter(qq_biological_sex == "Female") %>% pull(!!sym(var))
  male   <- data %>% filter(qq_biological_sex == "Male")   %>% pull(!!sym(var))
  f_agg <- mean_sem(female)
  m_agg <- mean_sem(male)
  female_str <- format_mean_sem(f_agg$mean, f_agg$sem)
  male_str   <- format_mean_sem(m_agg$mean, m_agg$sem)
  if (f_agg$n >= 2 && m_agg$n >= 2) {
    pval <- wilcox.test(female, male)$p.value
    if (is.na(pval) || !is.finite(pval)) {
      pval_str <- "N/A"
    } else if (round(pval, 4) == 1) {
      pval_str <- "1"
    } else {
      pval_str <- format(round(pval, 4), nsmall = 4)
    }
  } else {
    pval_str <- "N/A"
  }
  list(
    female = female_str,
    male = male_str,
    pvalue = pval_str,
    n_female = f_agg$n,
    n_male = m_agg$n
  )
}

# Build table body: one row per outcome, then 4 blocks of (Female, Male, p-value)
table_rows <- lapply(names(outcome_vars), function(v) {
  row <- c(Measure = outcome_vars[v])
  for (grp in group_levels) {
    block_data <- Table_4_Data %>% filter(qq_group == grp)
    bl <- build_block(block_data, v)
    row <- c(row,
             setNames(bl$female, paste0(grp, "_Female")),
             setNames(bl$male,   paste0(grp, "_Male")),
             setNames(bl$pvalue, paste0(grp, "_pvalue")))
  }
  row
})

# Collect n by group and sex for headers
n_headers <- Table_4_Data %>%
  count(qq_group, qq_biological_sex) %>%
  tidyr::pivot_wider(names_from = qq_biological_sex, values_from = n, values_fill = 0)

# Build matrix: rows = outcomes, columns = Measure, then for each group Female, Male, p-value
mat <- do.call(rbind, table_rows)
col_names <- c("Measure",
               paste0(rep(group_levels, each = 3), c("_Female", "_Male", "_pvalue")))
colnames(mat) <- col_names

# Simplify column names for display (we'll set proper headers below)
df_table <- as.data.frame(mat, stringsAsFactors = FALSE)

# Column structure: Measure, then Block1: Female (n=11), Male (n=3), p-value; Block2: ...; etc.
display_cols <- c("Measure")
for (grp in group_levels) {
  nf <- n_headers %>% filter(qq_group == grp) %>% pull(Female)
  nm <- n_headers %>% filter(qq_group == grp) %>% pull(Male)
  if (length(nf) == 0) nf <- 0
  if (length(nm) == 0) nm <- 0
  display_cols <- c(display_cols,
                   paste0(grp, "_Female"),
                   paste0(grp, "_Male"),
                   paste0(grp, "_pvalue"))
}

# Headers for flextable: one header row with group names spanning 3 columns each
ft <- flextable(df_table)

# Set header labels: Measure, then for each group "Female (n=xx)", "Male (n=xx)", "p-value"
header_labels <- c("Measure" = "Measure")
for (grp in group_levels) {
  nf <- n_headers %>% filter(qq_group == grp) %>% pull(Female)
  nm <- n_headers %>% filter(qq_group == grp) %>% pull(Male)
  if (length(nf) == 0) nf <- 0
  if (length(nm) == 0) nm <- 0
  header_labels <- c(header_labels,
                    setNames(paste0("Female (n=", nf, ")"), paste0(grp, "_Female")),
                    setNames(paste0("Male (n=", nm, ")"), paste0(grp, "_Male")),
                    setNames("p-value", paste0(grp, "_pvalue")))
}
ft <- set_header_labels(ft, values = header_labels)

# Add title and footnote
ft <- add_header_lines(ft, values = "Supplementary Table 4. Clinical and Biomarker Outcomes by COVID-19, mTBI, and Biological Sex.", top = TRUE)
ft <- add_footer_lines(ft, values = "p-values from Wilcoxon rank sum tests (Female vs Male within each group). Data are mean (SEM).")

# Theme and layout
ft <- theme_booktabs(ft)
ft <- width(ft, j = 1, width = 2.5)
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = 1, align = "left", part = "body")
ft <- align(ft, j = 1, align = "left", part = "header")

Supplementary_Table_4_Final <- ft
output_dir <- "supp_table_4/supp_table_4_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
save_as_docx(Supplementary_Table_4_Final, path = "supp_table_4/supp_table_4.docx")
write.csv(df_table, file.path(output_dir, "supp_table_4_export.csv"), row.names = FALSE)
