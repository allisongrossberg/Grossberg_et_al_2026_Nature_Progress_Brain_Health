# =============================================================================
# Supplementary Table 1 – Demographic and clinical characteristics
# =============================================================================
#
# Description: Summary table of demographics and clinical variables by study
#   group. Continuous: mean (SE); categorical: n (%). Tests: Wilcoxon/Kruskal–
#   Wallis for continuous; chi-square or Fisher for categorical. Exported as
#   Word and CSV.
#
# Prerequisites: supp_table_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: supp_table_1/supp_table_1.docx; supp_table_1_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) %in% paste0("supp_table_", 1:5)) setwd("..")
input_path <- "supp_table_1/supp_table_1_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
if (!file.exists(input_path)) stop("Need COAST_Study_Data_Clean_Age_Groups_add_dates.csv in supp_table_1/supp_table_1_input_data/.")
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(input_path, stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(plotrix)
library(gtsummary)
library(flextable)
library(sjstats)
library(broom)

# Select variables and order as in manuscript Supplementary Table 1.
Supplementary_Table_1_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% dplyr::select(
  qq_group, recruitment_site, age_years, qq_biological_sex, bmi, qq_education, race_category,
  qq_tbi_num, covid_pos_test_num, qq_covid_number,
  average_years_since_tbi, average_years_since_covid,
  recode_sev_score_total_covid, recode_sev_score_total_tbi,
  recode_freq_score_total_covid, recode_freq_score_total_tbi,
  recode_overall_sev_score, recode_overall_freq_score,
  chronic_acute_overall, chronic_tbi, chronic_covid
)

Supplementary_Table_1_Data$qq_group <- factor(Supplementary_Table_1_Data$qq_group, levels = c(
  "COVID-19 (-) mTBI (-)",
  "COVID-19 (+) mTBI (-)",
  "COVID-19 (-) mTBI (+)",
  "COVID-19 (+) mTBI (+)"
)
)

Supplementary_Table_1_Data$covid_pos_test_num <- as.numeric(Supplementary_Table_1_Data$covid_pos_test_num)
Supplementary_Table_1_Data$average_years_since_tbi <- as.numeric(Supplementary_Table_1_Data$average_years_since_tbi)
Supplementary_Table_1_Data$average_years_since_covid <- as.numeric(Supplementary_Table_1_Data$average_years_since_covid)

output_dir <- "supp_table_1/supp_table_1_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: for each variable, run appropriate test and return test name, statistic, df, p-value.
correct_stat_test <- function(data, variable, by, ...) {
  if (is.numeric(data[[variable]])) {
    # For continuous variables, use Kruskal-Wallis or Wilcoxon
    groups <- unique(data[[by]])
    if (length(groups) == 2) {
      # Wilcoxon test for two groups
      test_result <- wilcox.test(data[[variable]] ~ data[[by]])
      test_name <- "Wilcox"
      stat <- test_result$statistic
      df <- length(data[[variable]]) - 2  # n - 2 for Wilcoxon
    } else {
      # Kruskal-Wallis test for more than two groups
      test_result <- kruskal.test(data[[variable]] ~ data[[by]])
      test_name <- "Kruskal"
      stat <- test_result$statistic
      df <- test_result$parameter  # Already provided for Kruskal-Wallis
    }
    p_value <- signif(test_result$p.value, 5)  # Round to 5 significant figures
  } else {
    # For categorical variables, use Chi-squared or Fisher's exact test
    cont_table <- table(data[[by]], data[[variable]])
    expected <- chisq.test(cont_table)$expected
    if (any(expected < 5)) {
      # Fisher's exact test if any expected cell count is less than 5
      test_result <- fisher.test(cont_table, simulate.p.value = TRUE)
      test_name <- "Fisher's"
      stat <- NA
      df <- NA  # Fisher's exact test doesn't have df
    } else {
      # Chi-squared test otherwise
      test_result <- chisq.test(cont_table)
      test_name <- "χ²"
      stat <- test_result$statistic
      df <- test_result$parameter
    }
    p_value <- signif(test_result$p.value, 5)  # Round to 5 significant figures
  }
  
  # Return a tibble with the results
  tibble(
    test_name = test_name,
    statistic = unname(stat),
    df = unname(df),
    p_value = p_value
  )
}

# Build summary table (gtsummary) with custom test column.
Supplementary_Table_1 <- Supplementary_Table_1_Data %>%
  tbl_summary(
    by = qq_group,
    statistic = list(all_continuous() ~ "{mean} ({std.error})", all_categorical() ~ "{n} ({p}%)"),
    missing = "no",
    digits = list(all_categorical() ~ c(0, 1), all_continuous() ~ c(1, 1, 1)),
    label = list(recruitment_site ~ "Recruitment Site",
                 age_years ~ "Age (Years)",
                 qq_biological_sex ~ "Biological Sex",
                 bmi ~ "Body Mass Index (BMI)",
                 qq_education ~ "Education Level",
                 race_category ~ "Race",
                 qq_tbi_num ~ "Average Number of mTBIs",
                 covid_pos_test_num ~ "Number of Positive COVID-19 Incidences",
                 qq_covid_number ~ "Number of Suspected COVID-19 Incidences",
                 average_years_since_tbi ~ "Average Years Since Last mTBI",
                 average_years_since_covid ~ "Average Years Since Last COVID-19 Incidence",
                 recode_sev_score_total_covid ~ "Total COVID-19 Symptom Severity (Across All Incidences)",
                 recode_sev_score_total_tbi ~ "Total mTBI Symptom Severity (Across All Injuries)",
                 recode_freq_score_total_covid ~ "Total COVID-19 Symptom Frequency (Across All Incidences)",
                 recode_freq_score_total_tbi ~ "Total mTBI Symptom Frequency (Across All Injuries)",
                 recode_overall_sev_score ~ "Total Symptom Severity Score",
                 recode_overall_freq_score ~ "Total Frequency Severity Score",
                 chronic_acute_overall ~ "Chronic vs. Acute Symptoms",
                 chronic_tbi ~ "Chronic vs. Acute mTBI Symptoms",
                 chronic_covid ~ "Chronic vs. Acute COVID-19 Symptoms")) %>%
  add_stat(fns = everything() ~ correct_stat_test) %>%
  modify_header(
    test_name = "**Test**",
    statistic = "**Test Statistic**",
    df = "**Degrees of Freedom**",
    p_value = "**p-value**"
  ) %>%
  bold_labels() %>%
  modify_footnote(
    update = starts_with("stat_") ~
      "Mean (SE) for continuous variables; n (%) for categorical variables"
  )

Supplementary_Table_1_Final <- as_flex_table(Supplementary_Table_1)
Supplementary_Table_1_Final

save_as_docx(Supplementary_Table_1_Final, path = "supp_table_1/supp_table_1.docx")
write.csv(Supplementary_Table_1$table_body, file.path(output_dir, "supp_table_1_export.csv"), row.names = FALSE)

