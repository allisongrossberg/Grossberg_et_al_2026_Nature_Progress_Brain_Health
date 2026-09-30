# =============================================================================
# Supplementary Table 3 – Vaccination details and clinical/biomarker outcomes
#   by vaccination status, within each study group
# =============================================================================
#
# Description: A single merged table of vaccination dose count/vaccine type
#   and clinical/biomarker outcomes, both stratified by study group and then
#   by vaccination status (8 group x status columns throughout).
#
#   The "Received COVID-19 Vaccination" row is omitted since it's circular
#   with the column split (already shown in each column's "N ="). Dose
#   count and vaccine type are undefined for unvaccinated participants, so
#   their columns correctly show 0 (0%) or an empty cell.
#
#   No significance testing: with 0-1 unvaccinated participants in 3 of 4
#   groups, a formal test would not be meaningful.
#
# Prerequisites: supp_table_3_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv;
#   supp_table_3_input_data/vaccination_full.csv (or figure_1_input_data/
#   COAST_Data_Y1_8_1_24.csv + COAST_Headers_Y1_8_1_24.csv, to regenerate it);
#   figure_1 outputs: figure_1f (NeuroQOL), figure_1i (N3PA), figure_1h (PHQ8),
#   figure_1e (WAI) in figure_1/figure_1_output_data/.
#
# Output: supp_table_3/supp_table_3.docx; supp_table_3_output_data/*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) %in% paste0("supp_table_", 1:5)) setwd("..")
input_dir  <- "supp_table_3/supp_table_3_input_data"
output_dir <- "supp_table_3/supp_table_3_output_data"
docx_dir   <- "supp_table_3"
fig1_out   <- "figure_1/figure_1_output_data"

library(dplyr)
library(gtsummary)
library(flextable)
library(plotrix)  # provides std.error(), used by tbl_summary's "{std.error}" stat
library(readr)

path_main      <- file.path(input_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv")
path_vacc_full <- file.path(input_dir, "vaccination_full.csv")
path_nq   <- file.path(fig1_out, "figure_1f_NeuroQOL_plot_data.csv")
path_n3pa <- file.path(fig1_out, "figure_1i_N3PA_Abeta_plot_data.csv")
path_phq8 <- file.path(fig1_out, "figure_1h_PHQ8_plot_data.csv")
path_wai  <- file.path(fig1_out, "figure_1e_WAI_plot_data.csv")

# Vaccination fields are not in the cleaned pipeline output - extract once
# from the raw Qualtrics export and cache.
if (!file.exists(path_vacc_full)) {
  path_full_data    <- "figure_1/figure_1_input_data/COAST_Data_Y1_8_1_24.csv"
  path_full_headers <- "figure_1/figure_1_input_data/COAST_Headers_Y1_8_1_24.csv"
  if (!file.exists(path_full_data) || !file.exists(path_full_headers)) {
    stop("Need COAST_Data_Y1_8_1_24.csv and COAST_Headers_Y1_8_1_24.csv (raw Qualtrics export, figure_1/figure_1_input_data) to build vaccination_full.csv.")
  }
  headers_df <- read_csv(path_full_headers, show_col_types = FALSE, n_max = 0)
  raw_dat <- suppressMessages(read_csv(path_full_data, show_col_types = FALSE, col_types = cols(.default = "c"), name_repair = "minimal"))
  colnames(raw_dat) <- colnames(headers_df)
  main_ids <- read.csv(path_main, stringsAsFactors = FALSE) %>% dplyr::select(participant_id)
  vacc_cols <- c("participant_id", "qq_covid_vaccination_status",
                 "qq_covid_vaccination_doses___1", "qq_covid_vaccination_doses___2",
                 "qq_covid_vaccination_doses___3", "qq_covid_vaccination_doses___4",
                 "qq_covid_vaccine_type_dose_1", "qq_covid_vaccine_type_dose_2", "qq_covid_vaccine_type_dose_3")
  vacc <- raw_dat %>%
    dplyr::filter(participant_id %in% main_ids$participant_id) %>%
    dplyr::distinct(participant_id, .keep_all = TRUE) %>%
    dplyr::select(dplyr::all_of(vacc_cols))
  write.csv(vacc, path_vacc_full, row.names = FALSE)
}

missing <- c(
  if (!file.exists(path_main))      "COAST_Study_Data_Clean_Age_Groups_add_dates.csv (in supp_table_3_input_data)",
  if (!file.exists(path_nq))   "figure_1f_NeuroQOL_plot_data.csv (run figure_1/figure_1f.R)",
  if (!file.exists(path_n3pa)) "figure_1i_N3PA_Abeta_plot_data.csv (run figure_1/figure_1i.R)",
  if (!file.exists(path_phq8)) "figure_1h_PHQ8_plot_data.csv (run figure_1/figure_1h.R)",
  if (!file.exists(path_wai))  "figure_1e_WAI_plot_data.csv (run figure_1/figure_1e.R)"
)
missing <- missing[!sapply(missing, is.null)]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(path_main, stringsAsFactors = FALSE)
Vaccination_Full         <- read.csv(path_vacc_full, stringsAsFactors = FALSE)
NeuroQOL_Data_join       <- read.csv(path_nq, stringsAsFactors = FALSE)
N3PA_SRX_Data_Clean_join <- read.csv(path_n3pa, stringsAsFactors = FALSE)
PHQ8_plot_data           <- read.csv(path_phq8, stringsAsFactors = FALSE)
WAI_plot_data            <- read.csv(path_wai, stringsAsFactors = FALSE)

group_levels <- c("COVID-19 (+) mTBI (+)", "COVID-19 (+) mTBI (-)",
                   "COVID-19 (-) mTBI (+)", "COVID-19 (-) mTBI (-)")

# "Checked"/"Unchecked" -> "Yes"/"No" for the per-dose checkbox columns.
dose_cols <- grep("^qq_covid_vaccination_doses___", names(Vaccination_Full), value = TRUE)
transform_checked <- function(x) ifelse(x == "Checked", "Yes", ifelse(x == "Unchecked", "No", NA))
Vaccination_Full[dose_cols] <- lapply(Vaccination_Full[dose_cols], transform_checked)

Table_3_Data <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
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
            by = c("participant_id", "qq_group")) %>%
  left_join(Vaccination_Full %>%
              dplyr::select(participant_id, qq_covid_vaccination_status, dplyr::all_of(dose_cols),
                            qq_covid_vaccine_type_dose_1, qq_covid_vaccine_type_dose_2, qq_covid_vaccine_type_dose_3),
            by = "participant_id") %>%
  dplyr::filter(!is.na(qq_covid_vaccination_status)) %>%
  dplyr::mutate(
    qq_group = factor(qq_group, levels = group_levels),
    qq_covid_vaccination_status = factor(qq_covid_vaccination_status, levels = c("Yes", "No"),
                                          labels = c("Vaccinated", "Unvaccinated"))
  )

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(Table_3_Data, file.path(output_dir, "supp_table_3_input_data.csv"), row.names = FALSE)

continuous_vars <- c("recode_overall_neuro_psych_sev_score", "recode_overall_neuro_psych_freq_score",
                      "qq_phq8_average_score", "qq_eq5d_index_score", "qq_wai_2",
                      "Neuro_QOL_TScore", "Abeta42_Abeta40_Ratio_R", "vib_anti_abeta_1_42_igg_iga")

var_labels <- list(
  qq_covid_vaccination_doses___1 ~ "Received Vaccination Dose 1",
  qq_covid_vaccination_doses___2 ~ "Received Vaccination Dose 2",
  qq_covid_vaccination_doses___3 ~ "Received Vaccination Dose 3",
  qq_covid_vaccination_doses___4 ~ "Received Vaccination Dose 4",
  qq_covid_vaccine_type_dose_1 ~ "Vaccine Type, Dose 1",
  qq_covid_vaccine_type_dose_2 ~ "Vaccine Type, Dose 2",
  qq_covid_vaccine_type_dose_3 ~ "Vaccine Type, Dose 3",
  recode_overall_neuro_psych_sev_score ~ "Total Neuro/Psych Severity",
  recode_overall_neuro_psych_freq_score ~ "Total Neuro/Psych Frequency",
  qq_phq8_average_score ~ "PHQ-8",
  qq_eq5d_index_score ~ "EQ-5D",
  qq_wai_2 ~ "WAI-2",
  Neuro_QOL_TScore ~ "Neuro-QOL",
  Abeta42_Abeta40_Ratio_R ~ "Aβ42/Aβ40 Ratio",
  vib_anti_abeta_1_42_igg_iga ~ "Anti-Aβ1-42 IgG/IgA"
)

# Row order: vaccination details first, then clinical/biomarker outcomes.
# "Received COVID-19 Vaccination" itself is omitted (see header note - it is
# exactly the column split, already shown in each column's "N = " header).
row_order <- c(dose_cols, "qq_covid_vaccine_type_dose_1", "qq_covid_vaccine_type_dose_2",
               "qq_covid_vaccine_type_dose_3", continuous_vars)

Supplementary_Table_3 <- Table_3_Data %>%
  dplyr::select(-participant_id) %>%
  tbl_strata(
    strata = qq_group,
    .tbl_fun = ~ .x %>%
      tbl_summary(
        by = qq_covid_vaccination_status,
        include = dplyr::all_of(row_order),
        # Force continuous vs. categorical explicitly: gtsummary's
        # per-subgroup auto-detection otherwise miscategorizes low-unique-
        # value scores as categorical in the smallest strata.
        type = list(dplyr::all_of(continuous_vars) ~ "continuous",
                    dplyr::all_of(dose_cols) ~ "dichotomous",
                    dplyr::all_of(c("qq_covid_vaccine_type_dose_1",
                                     "qq_covid_vaccine_type_dose_2", "qq_covid_vaccine_type_dose_3")) ~ "categorical"),
        # Dose columns: show only the "Yes" level (one row per dose, % who
        # received it) rather than both Yes/No rows -- Yes/No are
        # complementary so showing both is redundant. "dichotomous" (not
        # "categorical") is required for the value= override below to apply.
        value = list(dplyr::all_of(dose_cols) ~ "Yes"),
        statistic = list(all_continuous() ~ "{mean} ({std.error})", all_categorical() ~ "{n} ({p}%)"),
        digits = list(all_continuous() ~ c(2, 2), all_categorical() ~ c(0, 1)),
        missing = "no",
        label = var_labels
      ),
    .header = "**{strata}**"
  ) %>%
  bold_labels() %>%
  modify_footnote(
    update = everything() ~
      paste("Mean (SEM) for clinical/biomarker measures; n (%) for vaccination dose/type.",
            "\"Received COVID-19 Vaccination\" is omitted here as it is exactly the Vaccinated/Unvaccinated",
            "column split (see each column's N); dose and vaccine-type rows are correctly 0 (0%) or empty",
            "under Unvaccinated, since neither is defined for someone who was never vaccinated.",
            "Descriptive statistics only; no significance testing performed",
            "(0-1 unvaccinated participants in 3 of 4 groups precludes a meaningful comparison).")
  )

Supplementary_Table_3_Final <- as_flex_table(Supplementary_Table_3)
Supplementary_Table_3_Final

save_as_docx(Supplementary_Table_3_Final, path = file.path(docx_dir, "supp_table_3.docx"))
write.csv(Supplementary_Table_3$table_body, file.path(output_dir, "supp_table_3_export.csv"), row.names = FALSE)
