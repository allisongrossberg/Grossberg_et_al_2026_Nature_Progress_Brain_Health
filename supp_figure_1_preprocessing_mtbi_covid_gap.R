# =============================================================================
# Supplementary Figure 1h/1i preprocessing – mTBI-to-COVID gap
# =============================================================================
#
# Description: For each double-exposed (COVID-19 (+) mTBI (+)) participant,
#   computes the gap in years between their most recent mTBI before their
#   first COVID-19 infection and that infection -- how long the injury had
#   to resolve before COVID-19 hit. This is anchored to the COVID-19
#   infection date, unlike supp_figure_1f's "years since exposure as of the
#   study visit." In this cohort, 53 of 54 double-exposed participants had
#   their mTBI before COVID-19 (only 2 had a later mTBI), so only the
#   mTBI-before-COVID-19 direction has enough data to analyze.
#
#   Recomputed from raw event dates (not retained in the cleaned pipeline
#   CSV). mTBI dates: month name + year (day imputed as the 15th). COVID-19
#   dates: month/day numeric, year coded as years-since-1900 (e.g. "122" = 2022).
#
# Prerequisites: COAST_Data_Y1_8_1_24.csv, COAST_Headers_Y1_8_1_24.csv (raw
#   Qualtrics export, project root); figure_1/figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv.
#
# Output: supp_figure_1/supp_figure_1_input_data/mtbi_covid_gap.csv
#   (participant_id, gap_years_mtbi_to_covid, n_tbi_before_covid,
#   n_tbi_after_covid)
# =============================================================================

if (basename(getwd()) == "supp_figure_1") setwd("..")

library(readr)
library(dplyr)

path_full_data <- "figure_1/figure_1_input_data/COAST_Data_Y1_8_1_24.csv"
path_full_headers <- "figure_1/figure_1_input_data/COAST_Headers_Y1_8_1_24.csv"
path_main <- "figure_1/figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
missing <- c(
  if (!file.exists(path_full_data)) path_full_data,
  if (!file.exists(path_full_headers)) path_full_headers,
  if (!file.exists(path_main)) path_main
)
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

headers_df <- read_csv(path_full_headers, show_col_types = FALSE, n_max = 0)
dat <- suppressMessages(read_csv(path_full_data, show_col_types = FALSE, col_types = cols(.default = "c"), name_repair = "minimal"))
colnames(dat) <- colnames(headers_df)

main <- read.csv(path_main, stringsAsFactors = FALSE) %>% dplyr::select(participant_id, qq_group)
dat <- dat %>% dplyr::filter(participant_id %in% main$participant_id) %>% dplyr::distinct(participant_id, .keep_all = TRUE)

month_to_num <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  ok <- is.na(n) & !is.na(x)
  n[ok] <- match(x[ok], month.name)
  n
}

# mTBI event dates: month name + year, day imputed as 15th.
n_tbi <- 20
tbi_dates <- matrix(as.Date(NA), nrow = nrow(dat), ncol = n_tbi)
for (i in seq_len(n_tbi)) {
  m <- month_to_num(dat[[paste0("qq_tbi_", i, "_month")]])
  y <- suppressWarnings(as.numeric(dat[[paste0("qq_tbi_", i, "_year")]]))
  tbi_dates[, i] <- as.Date(ifelse(is.na(m) | is.na(y), NA, paste(y, m, "15", sep = "-")), format = "%Y-%m-%d")
}
class(tbi_dates) <- "Date"

# COVID-19 positive-test event dates: numeric month/day; year coded as
# years-since-1900 (verified: sample raw values "121"/"122"/"123" ->
# 2021/2022/2023).
n_covid <- 10
covid_dates <- matrix(as.Date(NA), nrow = nrow(dat), ncol = n_covid)
for (i in seq_len(n_covid)) {
  m <- suppressWarnings(as.numeric(dat[[paste0("qq_covid_", i, "_pos_test_month")]]))
  d <- suppressWarnings(as.numeric(dat[[paste0("qq_covid_", i, "_pos_test_day")]]))
  y <- suppressWarnings(as.numeric(dat[[paste0("qq_covid_", i, "_pos_test_year")]])) + 1900
  covid_dates[, i] <- as.Date(ifelse(is.na(m) | is.na(d) | is.na(y), NA, paste(y, m, d, sep = "-")), format = "%Y-%m-%d")
}
class(covid_dates) <- "Date"

first_covid <- as.Date(apply(covid_dates, 1, function(r) { v <- r[!is.na(r)]; if (length(v) == 0) NA else min(v) }), origin = "1970-01-01")

# For each participant, the gap = first_covid - (most recent mTBI strictly
# before first_covid), in years. Undefined if no mTBI precedes first_covid.
gap_years <- rep(NA_real_, nrow(dat))
n_tbi_before <- rep(NA_integer_, nrow(dat))
n_tbi_after <- rep(NA_integer_, nrow(dat))
for (i in seq_len(nrow(dat))) {
  fc <- first_covid[i]
  row <- tbi_dates[i, ]
  if (is.na(fc) || all(is.na(row))) next
  before <- row[!is.na(row) & row < fc]
  after  <- row[!is.na(row) & row >= fc]
  n_tbi_before[i] <- length(before)
  n_tbi_after[i]  <- length(after)
  if (length(before) > 0) {
    gap_years[i] <- as.numeric(fc - max(before)) / 365.25
  }
}

out <- data.frame(
  participant_id = dat$participant_id,
  gap_years_mtbi_to_covid = gap_years,
  n_tbi_before_covid = n_tbi_before,
  n_tbi_after_covid = n_tbi_after
)

output_dir <- "supp_figure_1/supp_figure_1_input_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(output_dir, "mtbi_covid_gap.csv"), row.names = FALSE)

cat("Double-positive participants:", sum(dat$qq_group == "COVID-19 (+) mTBI (+)"), "\n")
cat("...with a defined mTBI-to-COVID gap (mTBI before first COVID):",
    sum(!is.na(out$gap_years_mtbi_to_covid) & dat$qq_group == "COVID-19 (+) mTBI (+)"), "\n")
cat("...with an mTBI on/after their first COVID (reverse-direction candidates, not separately testable - see header):",
    sum(!is.na(out$n_tbi_after_covid) & out$n_tbi_after_covid > 0 & dat$qq_group == "COVID-19 (+) mTBI (+)"), "\n")
cat("Saved:", file.path(output_dir, "mtbi_covid_gap.csv"), "\n")
