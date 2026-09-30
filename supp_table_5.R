# =============================================================================
# Supplementary Table 5A & 5B – Symptom severity (COVID vs mTBI)
# =============================================================================
#
# Description: Table 5A: COVID-19 symptom severity (COVID-19 (+) mTBI (-) vs
#   COVID-19 (+) mTBI (+)). Table 5B: mTBI symptom severity (COVID-19 (-) mTBI (+)
#   vs COVID-19 (+) mTBI (+)). Mean (SEM); Wilcoxon rank-sum. Exported as Word and CSV.
#
# Prerequisites: supp_table_5_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: supp_table_5/supp_table_5A.docx, supp_table_5/supp_table_5B.docx;
#   supp_table_5_output_data/supp_table_5A_export.csv, supp_table_5B_export.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) %in% paste0("supp_table_", 1:5)) setwd("..")
input_path <- "supp_table_5/supp_table_5_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
if (!file.exists(input_path)) stop("Need COAST_Study_Data_Clean_Age_Groups_add_dates.csv in supp_table_5/supp_table_5_input_data/.")
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(input_path, stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(tibble)
library(plotrix)   # std.error() for SEM in table cells
library(gtsummary)
library(flextable)

# Custom stat: Wilcoxon rank-sum, return test name, W statistic, df (empty), p-value (numeric for bold_p)
wilcoxon_stat <- function(data, variable, by, ...) {
  x <- data[[variable]]
  g <- data[[by]]
  ok <- !is.na(x) & !is.na(g)
  if (sum(ok) < 2) return(tibble(test_name = "Wilcoxon", statistic = NA_real_, df = NA_real_, p.value = NA_real_))
  res <- wilcox.test(x[ok] ~ g[ok], exact = FALSE)
  tibble(test_name = "Wilcoxon", statistic = unname(res$statistic), df = NA_real_, p.value = res$p.value)
}

output_dir <- "supp_table_5/supp_table_5_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Table 5A: COVID-19 (+) mTBI (-) vs COVID-19 (+) mTBI (+) — COVID symptom severity (Total)
# -----------------------------------------------------------------------------
Table_5A_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  filter(qq_group %in% c("COVID-19 (+) mTBI (-)", "COVID-19 (+) mTBI (+)")) %>%
  dplyr::select(qq_group,
                recode_sev_score_total_covid_fatigue,
                recode_sev_score_total_covid_headache,
                recode_sev_score_total_covid_insomnia_sleep_problems,
                recode_sev_score_total_covid_drowsiness,
                recode_sev_score_total_covid_photophobia_phonophobia_se,
                recode_sev_score_total_covid_brain_fog,
                recode_sev_score_total_covid_confusion,
                recode_sev_score_total_covid_memory_problems,
                recode_sev_score_total_covid_difficulty_concentrating,
                recode_sev_score_total_covid_difficulty_finding_words,
                recode_sev_score_total_covid_paresthesia,
                recode_sev_score_total_covid_los,
                recode_sev_score_total_covid_lot,
                recode_sev_score_total_covid_dizziness_lightheadedness,
                recode_sev_score_total_covid_difficulty_balancing)

Table_5A_Data$qq_group <- factor(Table_5A_Data$qq_group,
  levels = c("COVID-19 (+) mTBI (-)", "COVID-19 (+) mTBI (+)"))

Supplementary_Table_5A <- Table_5A_Data %>%
  tbl_summary(
    by = qq_group,
    statistic = list(all_continuous() ~ "{mean} ({std.error})"),
    missing = "no",
    type = list(
      recode_sev_score_total_covid_fatigue = "continuous",
      recode_sev_score_total_covid_headache = "continuous",
      recode_sev_score_total_covid_insomnia_sleep_problems = "continuous",
      recode_sev_score_total_covid_drowsiness = "continuous",
      recode_sev_score_total_covid_photophobia_phonophobia_se = "continuous",
      recode_sev_score_total_covid_brain_fog = "continuous",
      recode_sev_score_total_covid_confusion = "continuous",
      recode_sev_score_total_covid_memory_problems = "continuous",
      recode_sev_score_total_covid_difficulty_concentrating = "continuous",
      recode_sev_score_total_covid_difficulty_finding_words = "continuous",
      recode_sev_score_total_covid_paresthesia = "continuous",
      recode_sev_score_total_covid_los = "continuous",
      recode_sev_score_total_covid_lot = "continuous",
      recode_sev_score_total_covid_dizziness_lightheadedness = "continuous",
      recode_sev_score_total_covid_difficulty_balancing = "continuous"
    ),
    digits = list(all_continuous() ~ c(2, 2)),
    label = list(
      recode_sev_score_total_covid_fatigue ~ "Fatigue",
      recode_sev_score_total_covid_headache ~ "Headache",
      recode_sev_score_total_covid_insomnia_sleep_problems ~ "Insomnia/Sleep Problems",
      recode_sev_score_total_covid_drowsiness ~ "Drowsiness",
      recode_sev_score_total_covid_photophobia_phonophobia_se ~ "Sensitivity to light/noise",
      recode_sev_score_total_covid_brain_fog ~ "Brain fog",
      recode_sev_score_total_covid_confusion ~ "Confusion",
      recode_sev_score_total_covid_memory_problems ~ "Memory problems",
      recode_sev_score_total_covid_difficulty_concentrating ~ "Difficulty concentrating",
      recode_sev_score_total_covid_difficulty_finding_words ~ "Difficulty finding words",
      recode_sev_score_total_covid_paresthesia ~ "Tingling or prickling sensation",
      recode_sev_score_total_covid_los ~ "Loss of smell",
      recode_sev_score_total_covid_lot ~ "Loss of taste",
      recode_sev_score_total_covid_dizziness_lightheadedness ~ "Dizziness/lightheadedness",
      recode_sev_score_total_covid_difficulty_balancing ~ "Difficulty balancing"
    )
  ) %>%
  add_stat(fns = everything() ~ wilcoxon_stat) %>%
  modify_header(
    label = "** **",
    test_name = "**Test**",
    statistic = "**Test Statistic**",
    df = "**Degrees of Freedom**",
    p.value = "**p-value**"
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05) %>%
  modify_footnote(update = starts_with("stat_") ~ "Mean (SEM). Wilcoxon rank-sum test.")

Supplementary_Table_5A_Final <- as_flex_table(Supplementary_Table_5A)
save_as_docx(Supplementary_Table_5A_Final, path = "supp_table_5/supp_table_5A.docx")

# -----------------------------------------------------------------------------
# Table 5B: COVID-19 (-) mTBI (+) vs COVID-19 (+) mTBI (+) — mTBI symptom severity (Total)
# -----------------------------------------------------------------------------
Table_5B_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  filter(qq_group %in% c("COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)")) %>%
  dplyr::select(qq_group,
                recode_sev_score_total_tbi_headache,
                recode_sev_score_total_tbi_balance_problems,
                recode_sev_score_total_tbi_dizziness,
                recode_sev_score_total_tbi_lightheadedness,
                recode_sev_score_total_tbi_fatigue,
                recode_sev_score_total_tbi_trouble_falling_asleep,
                recode_sev_score_total_tbi_sleeping_more,
                recode_sev_score_total_tbi_drowsiness,
                recode_sev_score_total_tbi_light_sensitivity,
                recode_sev_score_total_tbi_noise_sensitivity,
                recode_sev_score_total_tbi_irritability,
                recode_sev_score_total_tbi_feeling_frustrated_impatient,
                recode_sev_score_total_tbi_taking_longer_to_think,
                recode_sev_score_total_tbi_restlessness,
                recode_sev_score_total_tbi_sadness,
                recode_sev_score_total_tbi_nervousness_anxiousness,
                recode_sev_score_total_tbi_feeling_more_emotional,
                recode_sev_score_total_tbi_feeling_slowed_down,
                recode_sev_score_total_tbi_in_a_fog,
                recode_sev_score_total_tbi_difficulty_concentrating,
                recode_sev_score_total_tbi_difficulty_remembering,
                recode_sev_score_total_tbi_blurred_vision)

Table_5B_Data$qq_group <- factor(Table_5B_Data$qq_group,
  levels = c("COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))

Supplementary_Table_5B <- Table_5B_Data %>%
  tbl_summary(
    by = qq_group,
    statistic = list(all_continuous() ~ "{mean} ({std.error})"),
    missing = "no",
    type = list(
      recode_sev_score_total_tbi_headache = "continuous",
      recode_sev_score_total_tbi_balance_problems = "continuous",
      recode_sev_score_total_tbi_dizziness = "continuous",
      recode_sev_score_total_tbi_lightheadedness = "continuous",
      recode_sev_score_total_tbi_fatigue = "continuous",
      recode_sev_score_total_tbi_trouble_falling_asleep = "continuous",
      recode_sev_score_total_tbi_sleeping_more = "continuous",
      recode_sev_score_total_tbi_drowsiness = "continuous",
      recode_sev_score_total_tbi_light_sensitivity = "continuous",
      recode_sev_score_total_tbi_noise_sensitivity = "continuous",
      recode_sev_score_total_tbi_irritability = "continuous",
      recode_sev_score_total_tbi_feeling_frustrated_impatient = "continuous",
      recode_sev_score_total_tbi_taking_longer_to_think = "continuous",
      recode_sev_score_total_tbi_restlessness = "continuous",
      recode_sev_score_total_tbi_sadness = "continuous",
      recode_sev_score_total_tbi_nervousness_anxiousness = "continuous",
      recode_sev_score_total_tbi_feeling_more_emotional = "continuous",
      recode_sev_score_total_tbi_feeling_slowed_down = "continuous",
      recode_sev_score_total_tbi_in_a_fog = "continuous",
      recode_sev_score_total_tbi_difficulty_concentrating = "continuous",
      recode_sev_score_total_tbi_difficulty_remembering = "continuous",
      recode_sev_score_total_tbi_blurred_vision = "continuous"
    ),
    digits = list(all_continuous() ~ c(2, 2)),
    label = list(
      recode_sev_score_total_tbi_headache ~ "Headache",
      recode_sev_score_total_tbi_balance_problems ~ "Balance Problems",
      recode_sev_score_total_tbi_dizziness ~ "Dizziness",
      recode_sev_score_total_tbi_lightheadedness ~ "Lightheadedness",
      recode_sev_score_total_tbi_fatigue ~ "Fatigue",
      recode_sev_score_total_tbi_trouble_falling_asleep ~ "Trouble Falling Asleep",
      recode_sev_score_total_tbi_sleeping_more ~ "Sleeping More Than Usual",
      recode_sev_score_total_tbi_drowsiness ~ "Drowsiness",
      recode_sev_score_total_tbi_light_sensitivity ~ "Light sensitivity",
      recode_sev_score_total_tbi_noise_sensitivity ~ "Noise sensitivity",
      recode_sev_score_total_tbi_irritability ~ "Irritability",
      recode_sev_score_total_tbi_feeling_frustrated_impatient ~ "Feeling frustrated or impatient",
      recode_sev_score_total_tbi_taking_longer_to_think ~ "Taking longer to think",
      recode_sev_score_total_tbi_restlessness ~ "Restlessness",
      recode_sev_score_total_tbi_sadness ~ "Sadness",
      recode_sev_score_total_tbi_nervousness_anxiousness ~ "Nervousness/anxiousness",
      recode_sev_score_total_tbi_feeling_more_emotional ~ "Feeling more emotional than usual",
      recode_sev_score_total_tbi_feeling_slowed_down ~ "Feeling slowed down",
      recode_sev_score_total_tbi_in_a_fog ~ "Brain fog",
      recode_sev_score_total_tbi_difficulty_concentrating ~ "Difficulty concentrating",
      recode_sev_score_total_tbi_difficulty_remembering ~ "Difficulty remembering",
      recode_sev_score_total_tbi_blurred_vision ~ "Blurred vision"
    )
  ) %>%
  add_stat(fns = everything() ~ wilcoxon_stat) %>%
  modify_header(
    label = "** **",
    test_name = "**Test**",
    statistic = "**Statistic**",
    df = "**Degrees of Freedom**",
    p.value = "**p-value**"
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05) %>%
  modify_footnote(update = starts_with("stat_") ~ "Mean (SEM). Wilcoxon rank-sum test.")

Supplementary_Table_5B_Final <- as_flex_table(Supplementary_Table_5B)
save_as_docx(Supplementary_Table_5B_Final, path = "supp_table_5/supp_table_5B.docx")

tb_a <- Supplementary_Table_5A$table_body
tb_a_export <- as.data.frame(lapply(tb_a, function(col) if (is.list(col)) sapply(col, paste, collapse = ", ") else col))
write.csv(tb_a_export, file.path(output_dir, "supp_table_5A_export.csv"), row.names = FALSE)
tb_b <- Supplementary_Table_5B$table_body
tb_b_export <- as.data.frame(lapply(tb_b, function(col) if (is.list(col)) sapply(col, paste, collapse = ", ") else col))
write.csv(tb_b_export, file.path(output_dir, "supp_table_5B_export.csv"), row.names = FALSE)
