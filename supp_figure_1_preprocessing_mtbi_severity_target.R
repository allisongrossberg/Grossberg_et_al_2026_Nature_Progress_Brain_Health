# =============================================================================
# Supplementary Figure 1h/1i preprocessing – mTBI severity of target incidence
# =============================================================================
#
# Description: Isolates the symptom severity/frequency of specifically the
#   most recent mTBI before each double-exposed participant's first COVID-19
#   infection -- the same single injury that
#   supp_figure_1_preprocessing_mtbi_covid_gap.R's gap variable and
#   supp_figure_1h_preprocessing_mtbi_count.R's count variable describe.
#   Uses the same severity/frequency scoring formula as
#   figure_1_preprocessing_step_2_symptom_data.R (None..Severe -> 0..5,
#   Never..Always -> 1..5, max(during,after) - before), but restricted to
#   only this one target incidence rather than summed across every mTBI a
#   participant ever had.
#
# Prerequisites: COAST_Data_Y1_8_1_24.csv, COAST_Headers_Y1_8_1_24.csv (raw
#   Qualtrics export, project root); figure_1/figure_1_output_data/COAST_Study_Data_F100_Clean.csv
#   (raw per-incidence severity/frequency strings, pre-recode - output of
#   figure_1_preprocessing_step_1_data_cleaning.R); figure_1/figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv.
#
# Output: supp_figure_1/supp_figure_1_input_data/mtbi_severity_target_incidence.csv
#   (participant_id, target_incidence, tbi_severity_target_incidence,
#   tbi_frequency_target_incidence)
# =============================================================================

if (basename(getwd()) == "supp_figure_1") setwd("..")

library(readr)
library(dplyr)

path_full_data    <- "figure_1/figure_1_input_data/COAST_Data_Y1_8_1_24.csv"
path_full_headers <- "figure_1/figure_1_input_data/COAST_Headers_Y1_8_1_24.csv"
path_main <- "figure_1/figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
path_f100 <- "figure_1/figure_1_output_data/COAST_Study_Data_F100_Clean.csv"
missing <- c(
  if (!file.exists(path_full_data))    path_full_data,
  if (!file.exists(path_full_headers)) path_full_headers,
  if (!file.exists(path_main)) path_main,
  if (!file.exists(path_f100)) path_f100
)
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = "; "))

# -----------------------------------------------------------------------
# Step 1: identify each participant's target incidence - the most recent
# mTBI before their first COVID-19 infection (identical date logic to
# supp_figure_1_gap_preprocessing.R, extended here to also record which
# incidence INDEX that is, not just the gap in years).
# -----------------------------------------------------------------------
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

n_tbi <- 20
tbi_dates <- matrix(as.Date(NA), nrow = nrow(dat), ncol = n_tbi)
for (i in seq_len(n_tbi)) {
  m <- month_to_num(dat[[paste0("qq_tbi_", i, "_month")]])
  y <- suppressWarnings(as.numeric(dat[[paste0("qq_tbi_", i, "_year")]]))
  tbi_dates[, i] <- as.Date(ifelse(is.na(m) | is.na(y), NA, paste(y, m, "15", sep = "-")), format = "%Y-%m-%d")
}
class(tbi_dates) <- "Date"

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

target_incidence <- rep(NA_integer_, nrow(dat))
for (i in seq_len(nrow(dat))) {
  fc <- first_covid[i]
  row <- tbi_dates[i, ]
  if (is.na(fc) || all(is.na(row))) next
  before_idx <- which(!is.na(row) & row < fc)
  if (length(before_idx) == 0) next
  target_incidence[i] <- before_idx[which.max(row[before_idx])]
}

target_df <- data.frame(participant_id = dat$participant_id, target_incidence = target_incidence)

# -----------------------------------------------------------------------
# Step 2: per-incidence TBI severity/frequency change score, computed
# directly from the raw pre-recode strings in COAST_Study_Data_F100_Clean.csv
# (same value mapping and same before-vs-during/after formula as
# figure_1_preprocessing_step_2_symptom_data.R's value_mapping()/
# compare_lists()/symptom_indicence_score(), applied here only to each
# participant's own target incidence rather than summed across all of them).
# -----------------------------------------------------------------------
f100 <- read.csv(path_f100, stringsAsFactors = FALSE)
f100 <- f100 %>% dplyr::left_join(target_df, by = "participant_id")

sev_map  <- c("None" = 0, "Mild" = 1, "Mild-Moderate" = 2, "Moderate" = 3, "Moderate-Severe" = 4, "Severe" = 5)
freq_map <- c("Never or almost never had/have symptom" = 1, "Sometimes had the symptom" = 2,
              "Often had the symptom" = 3, "Frequently had the symptom" = 4, "Always had the symptom" = 5)

# Same 25-symptom subset as tbi_neuro_psych_symptom_scores_list in
# figure_1_preprocessing_step_2_symptom_data.R (lines 578, 1036) - the
# "neuro_psych" composite deliberately excludes nausea, vomiting, and pain
# from the full 28-symptom TBI list (those are presumably part of a separate
# physical-symptom composite, not part of recode_tbi_neuro_psych_sev/freq_score).
tbi_symptoms <- c("headache", "balance_problems", "dizziness", "lightheadedness",
                   "fatigue", "trouble_falling_asleep", "sleeping_more", "sleeping_less", "drowsiness",
                   "light_sensitivity", "noise_sensitivity", "irritability", "feeling_frustrated_impatient",
                   "taking_longer_to_think", "restlessness", "sadness", "nervousness_anxiousness",
                   "feeling_more_emotional", "numbness_tingling", "feeling_slowed_down", "in_a_fog",
                   "difficulty_concentrating", "difficulty_remembering", "blurred_vision", "double_vision")

# Per-symptom score for ONE incidence: max(during, after) - before; NA if any
# of the three timepoints is missing (same convention as compare_lists() /
# symptom_indicence_score() in figure_1_preprocessing_step_2_symptom_data.R).
map_lookup <- function(map, x) {
  if (is.na(x)) return(NA_real_)
  val <- unname(map[as.character(x)])
  if (length(val) == 0 || is.na(val)) return(NA_real_)
  val
}

incidence_symptom_score <- function(row, i, symptom, kind, map) {
  col <- function(tp) paste0("qq_tbi_", i, "_", kind, "_", tp, "_", symptom)
  cols <- c(col("before"), col("during"), col("after"))
  if (!all(cols %in% names(row))) return(NA_real_)
  before <- map_lookup(map, row[[col("before")]])
  during <- map_lookup(map, row[[col("during")]])
  after  <- map_lookup(map, row[[col("after")]])
  if (is.na(before) || is.na(during) || is.na(after)) return(NA_real_)
  max(during, after) - before
}

# Total (across the 25 neuro_psych symptoms) for a participant's own target
# incidence. Matches recode_tbi_neuro_psych_sev/freq_score's own convention
# exactly (rowSums(..., na.rm = TRUE) with no all-NA guard at this final
# step - confirmed against it in the sanity check below): a participant with
# zero recorded data for every neuro_psych symptom scores 0 here, same as
# they would in the whole-life aggregate, not NA.
compute_target_total <- function(df, symptoms, kind, map) {
  vapply(seq_len(nrow(df)), function(idx) {
    i <- df$target_incidence[idx]
    if (is.na(i)) return(NA_real_)
    row <- df[idx, ]
    scores <- vapply(symptoms, function(s) incidence_symptom_score(row, i, s, kind, map), numeric(1))
    sum(scores, na.rm = TRUE)
  }, numeric(1))
}

f100$tbi_severity_target_incidence  <- compute_target_total(f100, tbi_symptoms, "sev",  sev_map)
f100$tbi_frequency_target_incidence <- compute_target_total(f100, tbi_symptoms, "freq", freq_map)

out <- f100 %>% dplyr::select(participant_id, target_incidence, tbi_severity_target_incidence, tbi_frequency_target_incidence)

output_dir <- "supp_figure_1/supp_figure_1_input_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(output_dir, "mtbi_severity_target_incidence.csv"), row.names = FALSE)

cat("Participants with a defined target incidence: n =", sum(!is.na(out$target_incidence)), "\n")
cat("...with a defined target-incidence severity score:", sum(!is.na(out$tbi_severity_target_incidence)), "\n")
cat("...with a defined target-incidence frequency score:", sum(!is.na(out$tbi_frequency_target_incidence)), "\n")

# -----------------------------------------------------------------------
# Sanity check: for participants whose ONLY mTBI (ever, before or after
# COVID-19) is the target incidence itself, the target-incidence score
# must exactly equal the whole-life aggregate (recode_tbi_neuro_psych_sev_score/
# freq_score), since in that case there is nothing else being summed over.
# -----------------------------------------------------------------------
check_df <- read.csv(path_main, stringsAsFactors = FALSE) %>%
  dplyr::select(participant_id, qq_tbi_num, recode_tbi_neuro_psych_sev_score, recode_tbi_neuro_psych_freq_score) %>%
  dplyr::left_join(out, by = "participant_id") %>%
  dplyr::filter(qq_tbi_num == 1, !is.na(target_incidence))
cat("\nSanity check (n =", nrow(check_df), "participants with exactly 1 lifetime mTBI):\n")
cat("Severity matches whole-life aggregate exactly:",
    sum(check_df$tbi_severity_target_incidence == check_df$recode_tbi_neuro_psych_sev_score, na.rm = TRUE),
    "of", sum(!is.na(check_df$tbi_severity_target_incidence)), "\n")
cat("Frequency matches whole-life aggregate exactly:",
    sum(check_df$tbi_frequency_target_incidence == check_df$recode_tbi_neuro_psych_freq_score, na.rm = TRUE),
    "of", sum(!is.na(check_df$tbi_frequency_target_incidence)), "\n")
