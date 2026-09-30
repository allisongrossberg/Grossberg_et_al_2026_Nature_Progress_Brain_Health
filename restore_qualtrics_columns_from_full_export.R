# =============================================================================
# Restore Qualtrics columns from full export
# =============================================================================
#
# Rebuilds figure_1_input_data/COAST_Data_Y1_8_1_24.csv and
# COAST_Headers_Y1_8_1_24.csv (the trimmed files the analysis scripts read,
# and the only Qualtrics-derived files meant to leave this machine) from the
# full raw Qualtrics export in Qualtrics_Data_Export_Full/. The full export
# contains direct identifiers (name, email, exact date of birth) and is
# never trimmed in place; this script computes age_years from date of birth
# once here, then keeps only the columns listed below -- date of birth,
# name, and email are never among them. Run this after adding a new column
# dependency to any script, then re-run figure_1's preprocessing steps.
#
# Usage: Rscript figure_1/restore_qualtrics_columns_from_full_export.R
#   (from the project root or from figure_1)
# =============================================================================

library(readr)
library(dplyr)
library(lubridate)

if (basename(getwd()) == "figure_1") setwd("..")
full_export <- "Qualtrics_Data_Export_Full"
input_dir   <- "figure_1/figure_1_input_data"

path_full_data    <- file.path(full_export, "COAST_Data_Y1_8_1_24.csv")
path_full_headers <- file.path(full_export, "COAST_Headers_Y1_8_1_24.csv")
if (!file.exists(path_full_data) || !file.exists(path_full_headers)) {
  stop("Full export files not found. Expected:\n  ", path_full_data, "\n  ", path_full_headers)
}
path_out_data    <- file.path(input_dir, "COAST_Data_Y1_8_1_24.csv")
path_out_headers <- file.path(input_dir, "COAST_Headers_Y1_8_1_24.csv")

message("Reading full export from: ", full_export)
dat <- read_csv(path_full_data, show_col_types = FALSE, col_types = cols(.default = "c"), name_repair = "minimal")
headers_df <- read_csv(path_full_headers, show_col_types = FALSE)
colnames(dat) <- colnames(headers_df)
n <- names(dat)

message("Full export has ", length(n), " columns.")

# Compute age_years from date of birth + survey start date here, once, so
# that date of birth itself never has to appear in the trimmed output
# (same formula figure_1_preprocessing_step_1_data_cleaning.R used to use
# directly on date of birth).
dob <- as.Date(paste(dat$qq_dob_month, dat$qq_dob_day, as.numeric(dat$qq_dob_year) + 1899, sep = "-"), format = "%m-%d-%Y")
start_date <- as.Date(sapply(strsplit(dat$qq_start_date, " "), `[`, 1), format = "%Y-%m-%d")
dat$age_years <- interval(as.POSIXct(dob), as.POSIXct(start_date)) %>% as.numeric('years')
n <- names(dat)

# Identifiers that must never appear in the trimmed output, regardless of
# which pattern below might otherwise match them.
identifier_cols <- c(
  "qq_dob_month", "qq_dob_day", "qq_dob_year",
  "cfha_consent_participant_dob", "cfha_hipaa_participant_dob",
  "first_name", "last_name", "qq_first_name", "qq_last_name",
  "cfha_consent_participant_first_name", "cfha_consent_participant_last_name",
  "cfha_hipaa_participant_first_name", "cfha_hipaa_participant_last_name",
  "sv_check_in_participant_first_name", "sv_check_in_participant_last_name",
  "new_illness_injury_first_name", "new_illness_injury_last_name",
  "qq_email_address", "cfha_consent_email_address", "cfha_consent_phone_number"
)

# Step 1 explicit columns
step1_explicit <- c(
  "participant_id", "qq_height", "age_years",
  "qq_start_date", "qq_race___1", "qq_race___2", "qq_race___3", "qq_race___4", "qq_race_other",
  "qq_weight_lbs", "vib_apoe_result",
  paste0("qq_pss_", 1:10),
  paste0("qq_neuroqol_", 1:28)
)
# Step 1 pattern-based (phq8, gad7, fss)
step1_pattern <- c(
  grep("_phq8_[0-9]+", n, value = TRUE),
  grep("_gad7_[0-9]+", n, value = TRUE),
  grep("_fss_[0-9]+", n, value = TRUE)
)
# Step 2 needs (pass-through from raw)
step2_need <- c(
  "qq_group", "qq_biological_sex", "recruitment_site", "qq_education", "qq_tbi_num", "qq_covid_number",
  paste0("qq_covid_", 1:10, "_test_results"),
  "qq_phq8_average_score", "qq_phq8_total_score", "qq_phq8_cat",
  "qq_eq5d_index_score", "qq_eq5d_mobility", "qq_eq5d_selfcare", "qq_eq5d_usual_activities",
  "qq_eq5d_pain_discomfort", "qq_eq5d_anxiety_depression",
  "qq_wai_1", "qq_wai_2", "qq_wai_3", "qq_wai_4", "qq_wai_5",
  "qq_fss_average_score", "qq_fss_total_score", "qq_fss_cat",
  "qq_gad7_average_score", "qq_gad7_total_score", "qq_gad7_cat",
  "vib_anti_abeta_1_42_igg_iga", "vib_anti_abeta_1_42_igm",
  # Needed by supp_table_3.R (vaccination status/dose count/vaccine type)
  "qq_covid_vaccination_status",
  paste0("qq_covid_vaccination_doses___", 1:4),
  paste0("qq_covid_vaccine_type_dose_", 1:3)
)
# All _sev_ and _freq_ columns (COVID and TBI symptom severity and frequency)
step2_sev  <- grep("_sev_", n, value = TRUE)
step2_freq <- grep("_freq_", n, value = TRUE)
# COVID pos_test date columns forced into month, day, year order so that
# step_2's unite(qq_covid_N_pos_test_month:qq_covid_N_pos_test_year) spans
# month-day-year and produces a string parseable by as.Date(..., "%m-%d-%Y").
step2_covid_date_order <- character(0)
for (i in 1:10) {
  for (part in c("month", "day", "year")) {
    nm <- paste0("qq_covid_", i, "_pos_test_", part)
    if (nm %in% n) step2_covid_date_order <- c(step2_covid_date_order, nm)
  }
}
# TBI date columns (month/year order is fine; step_2 handles TBI dates differently)
step2_tbi_dates <- grep("^qq_tbi_[0-9]+_(month|year)$", n, value = TRUE)
step2_dates <- c(step2_covid_date_order, step2_tbi_dates)
# COVID incidence: keep the full per-incident block, same as TBI below.
# A narrower COVID-specific pattern here would drop qq_covid_N_duration_*
# (chronic_covid's source columns) while keeping the equivalent TBI columns.
step2_extra <- n[grepl("^qq_covid_[0-9]+_", n) |
                  grepl("^qq_tbi_[0-9]+_", n)]
step2_extra <- step2_extra[!grepl("_sev_", step2_extra) & !grepl("_freq_", step2_extra)]

keep <- unique(c(step1_explicit, step1_pattern, step2_need, step2_sev, step2_freq, step2_dates, step2_extra))
keep <- keep[keep %in% n]
keep <- setdiff(keep, identifier_cols)  # defense in depth; see identifier_cols above

missing_need <- step2_need[!step2_need %in% keep]
if (length(missing_need) > 0) {
  stop("Full export is missing required columns: ", paste(missing_need, collapse = ", "))
}
if (length(step2_sev) == 0) {
  stop("Full export has no _sev_ columns; step 2 requires symptom severity data.")
}

dat_trimmed <- dat %>% select(all_of(keep))
stopifnot(!any(identifier_cols %in% names(dat_trimmed)))  # must never write an identifier column
headers_trimmed <- as.data.frame(matrix(names(dat_trimmed), nrow = 1))
names(headers_trimmed) <- names(dat_trimmed)

write.csv(dat_trimmed, path_out_data, row.names = FALSE)
write.csv(headers_trimmed, path_out_headers, row.names = FALSE)

message("Restored ", length(keep), " columns (", length(step2_sev), " _sev_, ", length(step2_freq), " _freq_).")
message("Written: ", path_out_data)
message("Written: ", path_out_headers)
message("Next: run figure_1_preprocessing_step_1_data_cleaning.R, then step_2.")
