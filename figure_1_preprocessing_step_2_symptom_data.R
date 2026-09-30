# =============================================================================
# COAST Study Year 1 – Symptom severity and derived variables (first 100 participants)
# =============================================================================
#
# Description: Builds on COAST_Study_Data_F100_Clean (from figure_1_preprocessing_step_1_data_cleaning.R)
#   to create COVID/TBI symptom severity and frequency scores, age groups, chronic/acute
#   classification, incidence-level dates, and the main analysis dataframe used by
#   etable scripts and Figure 1.
#
# Input:  COAST_Study_Data_F100_Clean (must be in environment; run cleaning script first)
# Output: COAST_Study_Data_Clean_Age_Groups_add_dates (main object for etables/figures)
#         COAST_Study_Data_Clean_only_covid_pos, COAST_Study_Data_Clean_only_covid_pos_check
#
# Usage:  Set working directory to figure_1. Run figure_1_preprocessing_step_1_data_cleaning.R
#         first, then this script (or load COAST_Study_Data_F100_Clean from figure_1_output_data).
#
# Packages: plyr (load before dplyr), zoo, gtools, lubridate, tidyr, dplyr, ggplot2,
#            ComplexHeatmap, circlize, tidyverse, dichromat, gridtext, tibble, stringr
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

# -----------------------------------------------------------------------------
# Packages (plyr before dplyr to avoid namespace conflicts)
# -----------------------------------------------------------------------------
library(plyr)
library(zoo)
library(gtools)
library(lubridate)
library(tidyr)
library(dplyr)
library(stringr)

# -----------------------------------------------------------------------------
# Paths (run from figure_1 or project root)
# -----------------------------------------------------------------------------
base_dir <- if (dir.exists("figure_1")) "figure_1" else "."

# -----------------------------------------------------------------------------
# Load cleaned data from script 1
# -----------------------------------------------------------------------------
COAST_Study_Data_F100_Clean <- read.csv(file.path(base_dir, "figure_1_output_data", "COAST_Study_Data_F100_Clean.csv"),
                                        stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Restrict COVID-19 data to positive-test incidences (set other incidence data to NA)
# -----------------------------------------------------------------------------
COAST_Study_Data_Clean_only_covid_pos <- COAST_Study_Data_F100_Clean %>%
  mutate_at(vars(starts_with("qq_covid_1")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_1_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_1_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_2")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_2_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_2_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_3")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_3_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_3_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_4")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_4_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_4_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_5")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_5_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_5_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_6")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_6_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_6_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_7")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_7_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_7_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_8")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_8_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_8_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_9")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_9_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_9_test_results), NA, .))) %>%
  mutate_at(vars(starts_with("qq_covid_10")), list(~ ifelse(COAST_Study_Data_F100_Clean$qq_covid_10_test_results == "No" | is.na(COAST_Study_Data_F100_Clean$qq_covid_10_test_results), NA, .)))

COAST_Study_Data_Clean_only_covid_pos_check <- COAST_Study_Data_Clean_only_covid_pos %>% select(qq_covid_1_test_results, qq_covid_2_test_results, qq_covid_3_test_results, qq_covid_4_test_results, qq_covid_5_test_results, 
                                                                                                qq_covid_6_test_results, qq_covid_7_test_results, qq_covid_8_test_results, qq_covid_9_test_results, qq_covid_10_test_results, qq_group, participant_id)

# Convert COVID test results from "Yes" to 1 
COAST_Study_Data_Clean_only_covid_pos_check <- COAST_Study_Data_Clean_only_covid_pos_check %>% 
  mutate_at(vars(starts_with("qq_covid_"),ends_with("_test_results")), list(~ ifelse(. == "Yes", 1, NA)))

# Sum number of positive COVID test results per participant 
COAST_Study_Data_Clean_only_covid_pos_check <- COAST_Study_Data_Clean_only_covid_pos_check %>% 
  mutate(sum_of_pos_tests = rowSums(select(., starts_with("qq_covid_"),ends_with("_test_results")), na.rm = TRUE))

table(COAST_Study_Data_Clean_only_covid_pos_check$sum_of_pos_tests)
table(COAST_Study_Data_Clean_only_covid_pos_check$sum_of_pos_tests[COAST_Study_Data_Clean_only_covid_pos_check$qq_group=="COVID-19 (+) mTBI (-)"])
table(COAST_Study_Data_Clean_only_covid_pos_check$sum_of_pos_tests[COAST_Study_Data_Clean_only_covid_pos_check$qq_group=="COVID-19 (+) mTBI (+)"])

# -----------------------------------------------------------------------------
# Recode COVID-19 and TBI symptom severity to numeric (_recode columns)
# -----------------------------------------------------------------------------
covid__tbi_symptom_sev_columns_to_recode <- grep("_sev_", names(COAST_Study_Data_Clean_only_covid_pos), value = TRUE)

for (col_name in covid__tbi_symptom_sev_columns_to_recode) {
  COAST_Study_Data_Clean_only_covid_pos[paste0("recode_", col_name)] <- ifelse(COAST_Study_Data_Clean_only_covid_pos[col_name] == "None", 0, 
                                                                               ifelse(COAST_Study_Data_Clean_only_covid_pos[col_name] == "Mild", 1, 
                                                                                      ifelse(COAST_Study_Data_Clean_only_covid_pos[col_name] == "Mild-Moderate", 2, 
                                                                                             ifelse(COAST_Study_Data_Clean_only_covid_pos[col_name] == "Moderate", 3, 
                                                                                                    ifelse(COAST_Study_Data_Clean_only_covid_pos[col_name] == "Moderate-Severe", 4, 
                                                                                                           ifelse(COAST_Study_Data_Clean_only_covid_pos[col_name] == "Severe", 5, 
                                                                                                                  ifelse(COAST_Study_Data_Clean_only_covid_pos[col_name] %in% c("None","Mild", "Mild-Moderate", "Moderate", "Moderate-Severe", "Severe"), COAST_Study_Data_Clean_only_covid_pos[col_name],NA)))))))
}

#####
#severity 
#functions
value_mapping <- function(list) {
  map <- c("None" = 0, "Mild" = 1, "Mild-Moderate" = 2, "Moderate" = 3, "Moderate-Severe" = 4, "Severe" = 5)
  new_list <- mapvalues(list, from = names(map), to = map)
  return(new_list)
}
#compare_lists func
compare_lists <- function(list1, list2) {
  # Initialize an empty result list
  result_list <- c()
  # Iterate through the lists
  for (i in 1:length(list1)) {
    #handle NAs
    if (is.na(list1[i]) | is.na(list2[i])) {
      result_list <- c(result_list, NA)
    } else if (list1[i] >= list2[i]) {
      # Append the higher value to the result list
      result_list <- c(result_list, list1[i])
    } else {
      # Append the higher value to the result list
      result_list <- c(result_list, list2[i])
    }
  }
  return(result_list)
}

#returns the scores with math applied
symptom_indicence_score <- function(list_of_lists, df) {
  for (list in list_of_lists) {
    before_list <- value_mapping(na.fill(df[list[1]], NA))
    during_list <- value_mapping(na.fill(df[list[2]], NA))
    after_list <- value_mapping(na.fill(df[list[3]], NA))
    larger_value_list <- compare_lists(during_list, after_list)
    final_list <- as.numeric(larger_value_list) - as.numeric(before_list)
    col_name <- list[4]
    df[col_name] <- as.numeric(final_list)
  }
  return(df)
}

#generate list of lists for symptom sev
generate_list_of_lists <- function(max_incidences, list) {
  final_list <- list()
  for (item in list) {
    for (i in seq(1, max_incidences)) {
      sub_list <- c()
      for (tp in c("before", "during", "after", "score")) {
        sub_list <- c(sub_list,sprintf(item, i, tp))
      }
      final_list <- c(final_list, list(sub_list))
    }
  }
  return(final_list)
}

#####

######

#symptom lists with placeholders -> list of lists
#All symptoms 

covid_all_symptom_list <- c("recode_qq_covid_%d_sev_%s_fatigue",
                            "recode_qq_covid_%d_sev_%s_insomnia_sleep_problems",
                            "recode_qq_covid_%d_sev_%s_mood_swings_or_irritability",
                            "recode_qq_covid_%d_sev_%s_drowsiness",
                            "recode_qq_covid_%d_sev_%s_reduced_blurred_vision",
                            "recode_qq_covid_%d_sev_%s_photophobia_phonophobia_se",
                            "recode_qq_covid_%d_sev_%s_brain_fog",
                            "recode_qq_covid_%d_sev_%s_confusion",
                            "recode_qq_covid_%d_sev_%s_memory_problems",
                            "recode_qq_covid_%d_sev_%s_difficulty_concentrating",
                            "recode_qq_covid_%d_sev_%s_delerium",
                            "recode_qq_covid_%d_sev_%s_difficulty_finding_words",
                            "recode_qq_covid_%d_sev_%s_paresthesia",
                            "recode_qq_covid_%d_sev_%s_headache",
                            "recode_qq_covid_%d_sev_%s_los",
                            "recode_qq_covid_%d_sev_%s_lot",
                            "recode_qq_covid_%d_sev_%s_dizziness_lightheadedness",
                            "recode_qq_covid_%d_sev_%s_difficulty_balancing",
                            "recode_qq_covid_%d_sev_%s_tremors",
                            "recode_qq_covid_%d_sev_%s_stroke",
                            "recode_qq_covid_%d_sev_%s_seizures",
                            "recode_qq_covid_%d_sev_%s_hypoacusis",
                            "recode_qq_covid_%d_sev_%s_numbness_hands_feet",
                            "recode_qq_covid_%d_sev_%s_hypoesthesia",  
                            "recode_qq_covid_%d_sev_%s_nasal_congestion", 
                            "recode_qq_covid_%d_sev_%s_sore_throat", 
                            "recode_qq_covid_%d_sev_%s_runny_nose", 
                            "recode_qq_covid_%d_sev_%s_ear_pain",
                            "recode_qq_covid_%d_sev_%s_cough", 
                            "recode_qq_covid_%d_sev_%s_sputum_production", 
                            "recode_qq_covid_%d_sev_%s_difficulty_breathing_sob", 
                            "recode_qq_covid_%d_sev_%s_hoarse_voice", 
                            "recode_qq_covid_%d_sev_%s_chest_pain_tightness", 
                            "recode_qq_covid_%d_sev_%s_chills", 
                            "recode_qq_covid_%d_sev_%s_swollen_lymph_nodes", 
                            "recode_qq_covid_%d_sev_%s_skipping_meals_appetite_loss", 
                            "recode_qq_covid_%d_sev_%s_sensitivity_heat_cold", 
                            "recode_qq_covid_%d_sev_%s_sweats", 
                            "recode_qq_covid_%d_sev_%s_white_red_purple_swollen_fingers_toes", 
                            "recode_qq_covid_%d_sev_%s_fever_feverish", 
                            "recode_qq_covid_%d_sev_%s_weight_loss", 
                            "recode_qq_covid_%d_sev_%s_tachycardia_arrhythmia_palpitations", 
                            "recode_qq_covid_%d_sev_%s_eye_soreness_discomfort", 
                            "recode_qq_covid_%d_sev_%s_abdominal_pain_stomachache", 
                            "recode_qq_covid_%d_sev_%s_diarrhea", 
                            "recode_qq_covid_%d_sev_%s_nausea_vomiting", 
                            "recode_qq_covid_%d_sev_%s_muscle_weakness", 
                            "recode_qq_covid_%d_sev_%s_muscle_pain_aches", 
                            "recode_qq_covid_%d_sev_%s_bone_and_joint_pain", 
                            "recode_qq_covid_%d_sev_%s_neck_back_pain")

tbi_all_symptom_list <- c("recode_qq_tbi_%d_sev_%s_headache",
                          "recode_qq_tbi_%d_sev_%s_nausea",
                          "recode_qq_tbi_%d_sev_%s_vomiting",
                          "recode_qq_tbi_%d_sev_%s_balance_problems",
                          "recode_qq_tbi_%d_sev_%s_dizziness",
                          "recode_qq_tbi_%d_sev_%s_lightheadedness",
                          "recode_qq_tbi_%d_sev_%s_fatigue",
                          "recode_qq_tbi_%d_sev_%s_trouble_falling_asleep",
                          "recode_qq_tbi_%d_sev_%s_sleeping_more",
                          "recode_qq_tbi_%d_sev_%s_sleeping_less",
                          "recode_qq_tbi_%d_sev_%s_drowsiness",
                          "recode_qq_tbi_%d_sev_%s_light_sensitivity",
                          "recode_qq_tbi_%d_sev_%s_noise_sensitivity",
                          "recode_qq_tbi_%d_sev_%s_irritability",
                          "recode_qq_tbi_%d_sev_%s_feeling_frustrated_impatient",
                          "recode_qq_tbi_%d_sev_%s_taking_longer_to_think",
                          "recode_qq_tbi_%d_sev_%s_restlessness",
                          "recode_qq_tbi_%d_sev_%s_sadness",
                          "recode_qq_tbi_%d_sev_%s_nervousness_anxiousness",
                          "recode_qq_tbi_%d_sev_%s_feeling_more_emotional",
                          "recode_qq_tbi_%d_sev_%s_numbness_tingling",
                          "recode_qq_tbi_%d_sev_%s_feeling_slowed_down",
                          "recode_qq_tbi_%d_sev_%s_in_a_fog",
                          "recode_qq_tbi_%d_sev_%s_difficulty_concentrating",
                          "recode_qq_tbi_%d_sev_%s_difficulty_remembering",
                          "recode_qq_tbi_%d_sev_%s_blurred_vision",
                          "recode_qq_tbi_%d_sev_%s_double_vision", 
                          "recode_qq_tbi_%d_sev_%s_pain")


#####
#generate lists of lists
covid_list_of_lists_all <- generate_list_of_lists(10, covid_all_symptom_list)
tbi_list_of_lists_all <- generate_list_of_lists(20, tbi_all_symptom_list)

#assign new df for testing
test_df_all <- COAST_Study_Data_Clean_only_covid_pos

#create single column for each symptom, each incidence
test_df_covid_all <- symptom_indicence_score(covid_list_of_lists_all, test_df_all)
test_df_covid_tbi_all <- symptom_indicence_score(tbi_list_of_lists_all, test_df_covid_all)


#aggregate across symptom
#this works for each symptom, make sure you switch tbi for covid when doing those symptoms and stick to the same naming convention. you will need to change the data frame that plugs nto the next step as well
# "Do this for every symptom for COVID and TBI"
test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_all %>%
  mutate(recode_sev_score_total_covid_fatigue = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fatigue")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fatigue"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fatigue")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_headache = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_headache")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_headache"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_headache")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_nasal_congestion = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nasal_congestion")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nasal_congestion"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nasal_congestion")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_sore_throat = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sore_throat")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sore_throat"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sore_throat")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_runny_nose = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_runny_nose")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_runny_nose"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_runny_nose")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_ear_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_ear_pain")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_ear_pain"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_ear_pain")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_cough = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_cough")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_cough"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_cough")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_sputum_production = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sputum_production")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sputum_production"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sputum_production")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_difficulty_breathing_sob = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_breathing_sob")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_breathing_sob"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_breathing_sob")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_hoarse_voice = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hoarse_voice")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hoarse_voice"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hoarse_voice")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_chest_pain_tightness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chest_pain_tightness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chest_pain_tightness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chest_pain_tightness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_chills = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chills")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chills"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chills")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_swollen_lymph_nodes = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_swollen_lymph_nodes")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_swollen_lymph_nodes"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_swollen_lymph_nodes")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_skipping_meals_appetite_loss = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_skipping_meals_appetite_loss")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_skipping_meals_appetite_loss"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_skipping_meals_appetite_loss")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_insomnia_sleep_problems = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_insomnia_sleep_problems")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_insomnia_sleep_problems"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_insomnia_sleep_problems")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_sensitivity_heat_cold = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sensitivity_heat_cold")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sensitivity_heat_cold"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sensitivity_heat_cold")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_sweats = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sweats")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sweats"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sweats")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_white_red_purple_swollen_fingers_toes = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_white_red_purple_swollen_fingers_toes")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_white_red_purple_swollen_fingers_toes"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_white_red_purple_swollen_fingers_toes")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_fever_feverish = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fever_feverish")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fever_feverish"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fever_feverish")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_mood_swings_irritability = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_mood_swings_irritability")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_mood_swings_irritability"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_mood_swings_irritability")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_weight_loss = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_weight_loss")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_weight_loss"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_weight_loss")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_drowsiness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_drowsiness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_drowsiness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_drowsiness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_tachycardia_arrhythmia_palpitations = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tachycardia_arrhythmia_palpitations")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tachycardia_arrhythmia_palpitations"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tachycardia_arrhythmia_palpitations")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_eye_soreness_discomfort = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_eye_soreness_discomfort")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_eye_soreness_discomfort"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_eye_soreness_discomfort")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_reduced_blurred_vision = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_reduced_blurred_vision")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_reduced_blurred_vision"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_reduced_blurred_vision")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_photophobia_phonophobia_se = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_photophobia_phonophobia_se")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_photophobia_phonophobia_se"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_photophobia_phonophobia_se")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_brain_fog = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_brain_fog")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_brain_fog"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_brain_fog")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_confusion = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_confusion")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_confusion"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_confusion")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_memory_problems = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_memory_problems")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_memory_problems"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_memory_problems")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_difficulty_concentrating = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_concentrating")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_concentrating"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_concentrating")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_delirium = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_delirium")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_delirium"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_delirium")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_difficulty_finding_words = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_finding_words")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_finding_words"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_finding_words")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_paresthesia = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_paresthesia")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_paresthesia"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_paresthesia")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_los = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_los")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_los"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_los")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_lot = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_lot")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_lot"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_lot")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_dizziness_lightheadedness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_dizziness_lightheadedness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_dizziness_lightheadedness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_dizziness_lightheadedness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_difficulty_balancing = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_balancing")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_balancing"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_balancing")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_tremors = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tremors")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tremors"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tremors")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_stroke = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_stroke")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_stroke"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_stroke")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_seizures = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_seizures")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_seizures"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_seizures")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_hypoacusis = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoacusis")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoacusis"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoacusis")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_numbness_hands_feet = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_numbness_hands_feet")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_numbness_hands_feet"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_numbness_hands_feet")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_hypoethesia = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoethesia")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoethesia"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoethesia")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_abdominal_pain_stomachache = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_abdominal_pain_stomachache")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_abdominal_pain_stomachache"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_abdominal_pain_stomachache")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_diarrhea = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_diarrhea")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_diarrhea"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_diarrhea")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_nausea_vomiting = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nausea_vomiting")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nausea_vomiting"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nausea_vomiting")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_muscle_weakness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_weakness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_weakness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_weakness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_muscle_pain_aches = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_pain_aches")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_pain_aches"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_pain_aches")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_bone_and_joint_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_bone_and_joint_pain")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_bone_and_joint_pain"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_bone_and_joint_pain")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_covid_neck_back_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_neck_back_pain")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_neck_back_pain"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_neck_back_pain")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_headache = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_headache")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_headache"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_headache")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_nausea = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nausea")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nausea"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nausea")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_vomiting = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_vomiting")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_vomiting"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_vomiting")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_balance_problems = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_balance_problems")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_balance_problems"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_balance_problems")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_dizziness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_dizziness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_dizziness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_dizziness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_lightheadedness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_lightheadedness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_lightheadedness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_lightheadedness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_fatigue = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_fatigue")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_fatigue"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_fatigue")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_trouble_falling_asleep = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_trouble_falling_asleep")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_trouble_falling_asleep"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_trouble_falling_asleep")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_sleeping_more = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_more")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_more"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_more")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_sleeping_less = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_less")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_less"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_less")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_drowsiness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_drowsiness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_drowsiness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_drowsiness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_light_sensitivity = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_light_sensitivity")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_light_sensitivity"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_light_sensitivity")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_noise_sensitivity = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_noise_sensitivity")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_noise_sensitivity"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_noise_sensitivity")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_irritability = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_irritability")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_irritability"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_irritability")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_feeling_frustrated_impatient = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_frustrated_impatient")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_frustrated_impatient"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_frustrated_impatient")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_taking_longer_to_think = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_taking_longer_to_think")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_taking_longer_to_think"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_taking_longer_to_think")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_restlessness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_restlessness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_restlessness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_restlessness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_sadness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sadness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sadness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sadness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_nervousness_anxiousness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nervousness_anxiousness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nervousness_anxiousness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nervousness_anxiousness")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_feeling_more_emotional = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_more_emotional")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_more_emotional"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_more_emotional")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_numbness_tingling = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_numbness_tingling")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_numbness_tingling"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_numbness_tingling")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_feeling_slowed_down = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_slowed_down")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_slowed_down"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_slowed_down")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_in_a_fog = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_in_a_fog")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_in_a_fog"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_in_a_fog")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_difficulty_concentrating = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_concentrating")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_concentrating"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_concentrating")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_difficulty_remembering = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_remembering")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_remembering"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_remembering")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_blurred_vision = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_blurred_vision")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_blurred_vision"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_blurred_vision")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_double_vision = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_double_vision")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_double_vision"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_double_vision")), na.rm=TRUE))) %>% 
  mutate(recode_sev_score_total_tbi_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_pain")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_pain"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_pain")), na.rm=TRUE))) 

#calculate count for every symptom
# "Do this for every symptom for COVID and TBI"
test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals %>%
  mutate(recode_covid_fatigue_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fatigue"))))) %>%
  mutate(recode_covid_headache_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_headache"))))) %>%
  mutate(recode_covid_nasal_congestion_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nasal_congestion"))))) %>%
  mutate(recode_covid_sore_throat_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sore_throat"))))) %>%
  mutate(recode_covid_runny_nose_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_runny_nose"))))) %>%
  mutate(recode_covid_ear_pain_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_ear_pain"))))) %>%
  mutate(recode_covid_cough_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_cough"))))) %>%
  mutate(recode_covid_sputum_production_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sputum_production"))))) %>%
  mutate(recode_covid_difficulty_breathing_sob_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_breathing_sob"))))) %>%
  mutate(recode_covid_hoarse_voice_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hoarse_voice"))))) %>%
  mutate(recode_covid_chest_pain_tightness_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chest_pain_tightness"))))) %>%
  mutate(recode_covid_chills_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_chills"))))) %>%
  mutate(recode_covid_swollen_lymph_nodes_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_swollen_lymph_nodes"))))) %>%
  mutate(recode_covid_skipping_meals_appetite_loss_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_skipping_meals_appetite_loss"))))) %>%
  mutate(recode_covid_insomnia_sleep_problems_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_insomnia_sleep_problems"))))) %>%
  mutate(recode_covid_sensitivity_heat_cold_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sensitivity_heat_cold"))))) %>%
  mutate(recode_covid_sweats_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_sweats"))))) %>%
  mutate(recode_covid_white_red_purple_swollen_fingers_toes_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_white_red_purple_swollen_fingers_toes"))))) %>%
  mutate(recode_covid_fever_feverish_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_fever_feverish"))))) %>%
  mutate(recode_covid_mood_swings_irritability_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_mood_swings_irritability"))))) %>%
  mutate(recode_covid_weight_loss_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_weight_loss"))))) %>%
  mutate(recode_covid_drowsiness_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_drowsiness"))))) %>%
  mutate(recode_covid_tachycardia_arrhythmia_palpitations_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tachycardia_arrhythmia_palpitations"))))) %>%
  mutate(recode_covid_eye_soreness_discomfort_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_eye_soreness_discomfort"))))) %>%
  mutate(recode_covid_reduced_blurred_vision_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_reduced_blurred_vision"))))) %>%
  mutate(recode_covid_photophobia_phonophobia_se_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_photophobia_phonophobia_se"))))) %>%
  mutate(recode_covid_brain_fog_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_brain_fog"))))) %>%
  mutate(recode_covid_confusion_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_confusion"))))) %>%
  mutate(recode_covid_memory_problems_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_memory_problems"))))) %>%
  mutate(recode_covid_difficulty_concentrating_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_concentrating"))))) %>%
  mutate(recode_covid_delirium_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_delirium"))))) %>%
  mutate(recode_covid_difficulty_finding_words_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_finding_words"))))) %>%
  mutate(recode_covid_paresthesia_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_paresthesia"))))) %>%
  mutate(recode_covid_los_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_los"))))) %>%
  mutate(recode_covid_lot_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_lot"))))) %>%
  mutate(recode_covid_dizziness_lightheadedness_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_dizziness_lightheadedness"))))) %>%
  mutate(recode_covid_difficulty_balancing_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_difficulty_balancing"))))) %>%
  mutate(recode_covid_tremors_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_tremors"))))) %>%
  mutate(recode_covid_stroke_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_stroke"))))) %>%
  mutate(recode_covid_seizures_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_seizures"))))) %>%
  mutate(recode_covid_hypoacusis_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoacusis"))))) %>%
  mutate(recode_covid_numbness_hands_feet_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_numbness_hands_feet"))))) %>%
  mutate(recode_covid_hypoethesia_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_hypoethesia"))))) %>%
  mutate(recode_covid_abdominal_pain_stomachache_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_abdominal_pain_stomachache"))))) %>%
  mutate(recode_covid_diarrhea_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_diarrhea"))))) %>%
  mutate(recode_covid_nausea_vomiting_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_nausea_vomiting"))))) %>%
  mutate(recode_covid_muscle_weakness_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_weakness"))))) %>%
  mutate(recode_covid_muscle_pain_aches_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_muscle_pain_aches"))))) %>%
  mutate(recode_covid_bone_and_joint_pain_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_bone_and_joint_pain"))))) %>%
  mutate(recode_covid_neck_back_pain_count = rowSums(!is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_neck_back_pain"))))) %>%
  mutate(recode_tbi_headache_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_headache"))))) %>%
  mutate(recode_tbi_nausea_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nausea"))))) %>%
  mutate(recode_tbi_vomiting_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_vomiting"))))) %>%
  mutate(recode_tbi_balance_problems_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_balance_problems"))))) %>%
  mutate(recode_tbi_dizziness_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_dizziness"))))) %>%
  mutate(recode_tbi_lightheadedness_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_lightheadedness"))))) %>%
  mutate(recode_tbi_fatigue_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_fatigue"))))) %>%
  mutate(recode_tbi_trouble_falling_asleep_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_trouble_falling_asleep"))))) %>%
  mutate(recode_tbi_sleeping_more_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_more"))))) %>%
  mutate(recode_tbi_sleeping_less_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sleeping_less"))))) %>%
  mutate(recode_tbi_drowsiness_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_drowsiness"))))) %>%
  mutate(recode_tbi_light_sensitivity_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_light_sensitivity"))))) %>%
  mutate(recode_tbi_noise_sensitivity_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_noise_sensitivity"))))) %>%
  mutate(recode_tbi_irritability_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_irritability"))))) %>%
  mutate(recode_tbi_feeling_frustrated_impatient_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_frustrated_impatient"))))) %>%
  mutate(recode_tbi_taking_longer_to_think_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_taking_longer_to_think"))))) %>%
  mutate(recode_tbi_restlessness_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_restlessness"))))) %>%
  mutate(recode_tbi_sadness_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_sadness"))))) %>%
  mutate(recode_tbi_nervousness_anxiousness_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_nervousness_anxiousness"))))) %>%
  mutate(recode_tbi_feeling_more_emotional_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_more_emotional"))))) %>%
  mutate(recode_tbi_numbness_tingling_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_numbness_tingling"))))) %>%
  mutate(recode_tbi_feeling_slowed_down_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_feeling_slowed_down"))))) %>%
  mutate(recode_tbi_in_a_fog_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_in_a_fog"))))) %>%
  mutate(recode_tbi_difficulty_concentrating_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_concentrating"))))) %>%
  mutate(recode_tbi_difficulty_remembering_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_difficulty_remembering"))))) %>%
  mutate(recode_tbi_blurred_vision_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_blurred_vision"))))) %>%
  mutate(recode_tbi_double_vision_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_double_vision"))))) %>%
  mutate(recode_tbi_pain_count = rowSums(!is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_pain")))))

#calculate average for every symptom
# "Do this for every symptom for COVID and TBI"
test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals %>%
  mutate(recode_covid_fatigue_sev_average = recode_sev_score_total_covid_fatigue / recode_covid_fatigue_count) %>%
  mutate(recode_covid_headache_sev_average = recode_sev_score_total_covid_headache / recode_covid_headache_count) %>%
  mutate(recode_covid_nasal_congestion_sev_average = recode_sev_score_total_covid_nasal_congestion / recode_covid_nasal_congestion_count) %>%
  mutate(recode_covid_sore_throat_sev_average = recode_sev_score_total_covid_sore_throat / recode_covid_sore_throat_count) %>%
  mutate(recode_covid_runny_nose_sev_average = recode_sev_score_total_covid_runny_nose / recode_covid_runny_nose_count) %>%
  mutate(recode_covid_ear_pain_sev_average = recode_sev_score_total_covid_ear_pain / recode_covid_ear_pain_count) %>%
  mutate(recode_covid_cough_sev_average = recode_sev_score_total_covid_cough / recode_covid_cough_count) %>%
  mutate(recode_covid_sputum_production_sev_average = recode_sev_score_total_covid_sputum_production / recode_covid_sputum_production_count) %>%
  mutate(recode_covid_difficulty_breathing_sob_sev_average = recode_sev_score_total_covid_difficulty_breathing_sob / recode_covid_difficulty_breathing_sob_count) %>%
  mutate(recode_covid_hoarse_voice_sev_average = recode_sev_score_total_covid_hoarse_voice / recode_covid_hoarse_voice_count) %>%
  mutate(recode_covid_chest_pain_tightness_sev_average = recode_sev_score_total_covid_chest_pain_tightness / recode_covid_chest_pain_tightness_count) %>%
  mutate(recode_covid_chills_sev_average = recode_sev_score_total_covid_chills / recode_covid_chills_count) %>%
  mutate(recode_covid_swollen_lymph_nodes_sev_average = recode_sev_score_total_covid_swollen_lymph_nodes / recode_covid_swollen_lymph_nodes_count) %>%
  mutate(recode_covid_skipping_meals_appetite_loss_sev_average = recode_sev_score_total_covid_skipping_meals_appetite_loss / recode_covid_skipping_meals_appetite_loss_count) %>%
  mutate(recode_covid_insomnia_sleep_problems_sev_average = recode_sev_score_total_covid_insomnia_sleep_problems / recode_covid_insomnia_sleep_problems_count) %>%
  mutate(recode_covid_sensitivity_heat_cold_sev_average = recode_sev_score_total_covid_sensitivity_heat_cold / recode_covid_sensitivity_heat_cold_count) %>%
  mutate(recode_covid_sweats_sev_average = recode_sev_score_total_covid_sweats / recode_covid_sweats_count) %>%
  mutate(recode_covid_white_red_purple_swollen_fingers_toes_sev_average = recode_sev_score_total_covid_white_red_purple_swollen_fingers_toes / recode_covid_white_red_purple_swollen_fingers_toes_count) %>%
  mutate(recode_covid_fever_feverish_sev_average = recode_sev_score_total_covid_fever_feverish / recode_covid_fever_feverish_count) %>%
  mutate(recode_covid_mood_swings_irritability_sev_average = recode_sev_score_total_covid_mood_swings_irritability / recode_covid_mood_swings_irritability_count) %>%
  mutate(recode_covid_weight_loss_sev_average = recode_sev_score_total_covid_weight_loss / recode_covid_weight_loss_count) %>%
  mutate(recode_covid_drowsiness_sev_average = recode_sev_score_total_covid_drowsiness / recode_covid_drowsiness_count) %>%
  mutate(recode_covid_tachycardia_arrhythmia_palpitations_sev_average = recode_sev_score_total_covid_tachycardia_arrhythmia_palpitations / recode_covid_tachycardia_arrhythmia_palpitations_count) %>%
  mutate(recode_covid_eye_soreness_discomfort_sev_average = recode_sev_score_total_covid_eye_soreness_discomfort / recode_covid_eye_soreness_discomfort_count) %>%
  mutate(recode_covid_reduced_blurred_vision_sev_average = recode_sev_score_total_covid_reduced_blurred_vision / recode_covid_reduced_blurred_vision_count) %>%
  mutate(recode_covid_photophobia_phonophobia_se_sev_average = recode_sev_score_total_covid_photophobia_phonophobia_se / recode_covid_photophobia_phonophobia_se_count) %>%
  mutate(recode_covid_brain_fog_sev_average = recode_sev_score_total_covid_brain_fog / recode_covid_brain_fog_count) %>%
  mutate(recode_covid_confusion_sev_average = recode_sev_score_total_covid_confusion / recode_covid_confusion_count) %>%
  mutate(recode_covid_memory_problems_sev_average = recode_sev_score_total_covid_memory_problems / recode_covid_memory_problems_count) %>%
  mutate(recode_covid_difficulty_concentrating_sev_average = recode_sev_score_total_covid_difficulty_concentrating / recode_covid_difficulty_concentrating_count) %>%
  mutate(recode_covid_delirium_sev_average = recode_sev_score_total_covid_delirium / recode_covid_delirium_count) %>%
  mutate(recode_covid_difficulty_finding_words_sev_average = recode_sev_score_total_covid_difficulty_finding_words / recode_covid_difficulty_finding_words_count) %>%
  mutate(recode_covid_paresthesia_sev_average = recode_sev_score_total_covid_paresthesia / recode_covid_paresthesia_count) %>%
  mutate(recode_covid_los_sev_average = recode_sev_score_total_covid_los / recode_covid_los_count) %>%
  mutate(recode_covid_lot_sev_average = recode_sev_score_total_covid_lot / recode_covid_lot_count) %>%
  mutate(recode_covid_dizziness_lightheadedness_sev_average = recode_sev_score_total_covid_dizziness_lightheadedness / recode_covid_dizziness_lightheadedness_count) %>%
  mutate(recode_covid_difficulty_balancing_sev_average = recode_sev_score_total_covid_difficulty_balancing / recode_covid_difficulty_balancing_count) %>%
  mutate(recode_covid_tremors_sev_average = recode_sev_score_total_covid_tremors / recode_covid_tremors_count) %>%
  mutate(recode_covid_stroke_sev_average = recode_sev_score_total_covid_stroke / recode_covid_stroke_count) %>%
  mutate(recode_covid_seizures_sev_average = recode_sev_score_total_covid_seizures / recode_covid_seizures_count) %>%
  mutate(recode_covid_hypoacusis_sev_average = recode_sev_score_total_covid_hypoacusis / recode_covid_hypoacusis_count) %>%
  mutate(recode_covid_numbness_hands_feet_sev_average = recode_sev_score_total_covid_numbness_hands_feet / recode_covid_numbness_hands_feet_count) %>%
  mutate(recode_covid_hypoethesia_sev_average = recode_sev_score_total_covid_hypoethesia / recode_covid_hypoethesia_count) %>%
  mutate(recode_covid_abdominal_pain_stomachache_sev_average = recode_sev_score_total_covid_abdominal_pain_stomachache / recode_covid_abdominal_pain_stomachache_count) %>%
  mutate(recode_covid_diarrhea_sev_average = recode_sev_score_total_covid_diarrhea / recode_covid_diarrhea_count) %>%
  mutate(recode_covid_nausea_vomiting_sev_average = recode_sev_score_total_covid_nausea_vomiting / recode_covid_nausea_vomiting_count) %>%
  mutate(recode_covid_muscle_weakness_sev_average = recode_sev_score_total_covid_muscle_weakness / recode_covid_muscle_weakness_count) %>%
  mutate(recode_covid_muscle_pain_aches_sev_average = recode_sev_score_total_covid_muscle_pain_aches / recode_covid_muscle_pain_aches_count) %>%
  mutate(recode_covid_bone_and_joint_pain_sev_average = recode_sev_score_total_covid_bone_and_joint_pain / recode_covid_bone_and_joint_pain_count) %>%
  mutate(recode_covid_neck_back_pain_sev_average = recode_sev_score_total_covid_neck_back_pain / recode_covid_neck_back_pain_count) %>%
  mutate(recode_tbi_headache_sev_average = recode_sev_score_total_tbi_headache / recode_tbi_headache_count) %>%
  mutate(recode_tbi_nausea_sev_average = recode_sev_score_total_tbi_nausea / recode_tbi_nausea_count) %>%
  mutate(recode_tbi_vomiting_sev_average = recode_sev_score_total_tbi_vomiting / recode_tbi_vomiting_count) %>%
  mutate(recode_tbi_balance_problems_sev_average = recode_sev_score_total_tbi_balance_problems / recode_tbi_balance_problems_count) %>%
  mutate(recode_tbi_dizziness_sev_average = recode_sev_score_total_tbi_dizziness / recode_tbi_dizziness_count) %>%
  mutate(recode_tbi_lightheadedness_sev_average = recode_sev_score_total_tbi_lightheadedness / recode_tbi_lightheadedness_count) %>%
  mutate(recode_tbi_fatigue_sev_average = recode_sev_score_total_tbi_fatigue / recode_tbi_fatigue_count) %>%
  mutate(recode_tbi_trouble_falling_asleep_sev_average = recode_sev_score_total_tbi_trouble_falling_asleep / recode_tbi_trouble_falling_asleep_count) %>%
  mutate(recode_tbi_sleeping_more_sev_average = recode_sev_score_total_tbi_sleeping_more / recode_tbi_sleeping_more_count) %>%
  mutate(recode_tbi_sleeping_less_sev_average = recode_sev_score_total_tbi_sleeping_less / recode_tbi_sleeping_less_count) %>%
  mutate(recode_tbi_drowsiness_sev_average = recode_sev_score_total_tbi_drowsiness / recode_tbi_drowsiness_count) %>%
  mutate(recode_tbi_light_sensitivity_sev_average = recode_sev_score_total_tbi_light_sensitivity / recode_tbi_light_sensitivity_count) %>%
  mutate(recode_tbi_noise_sensitivity_sev_average = recode_sev_score_total_tbi_noise_sensitivity / recode_tbi_noise_sensitivity_count) %>%
  mutate(recode_tbi_irritability_sev_average = recode_sev_score_total_tbi_irritability / recode_tbi_irritability_count) %>%
  mutate(recode_tbi_feeling_frustrated_impatient_sev_average = recode_sev_score_total_tbi_feeling_frustrated_impatient / recode_tbi_feeling_frustrated_impatient_count) %>%
  mutate(recode_tbi_taking_longer_to_think_sev_average = recode_sev_score_total_tbi_taking_longer_to_think / recode_tbi_taking_longer_to_think_count) %>%
  mutate(recode_tbi_restlessness_sev_average = recode_sev_score_total_tbi_restlessness / recode_tbi_restlessness_count) %>%
  mutate(recode_tbi_sadness_sev_average = recode_sev_score_total_tbi_sadness / recode_tbi_sadness_count) %>%
  mutate(recode_tbi_nervousness_anxiousness_sev_average = recode_sev_score_total_tbi_nervousness_anxiousness / recode_tbi_nervousness_anxiousness_count) %>%
  mutate(recode_tbi_feeling_more_emotional_sev_average = recode_sev_score_total_tbi_feeling_more_emotional / recode_tbi_feeling_more_emotional_count) %>%
  mutate(recode_tbi_numbness_tingling_sev_average = recode_sev_score_total_tbi_numbness_tingling / recode_tbi_numbness_tingling_count) %>%
  mutate(recode_tbi_feeling_slowed_down_sev_average = recode_sev_score_total_tbi_feeling_slowed_down / recode_tbi_feeling_slowed_down_count) %>%
  mutate(recode_tbi_in_a_fog_sev_average = recode_sev_score_total_tbi_in_a_fog / recode_tbi_in_a_fog_count) %>%
  mutate(recode_tbi_difficulty_concentrating_sev_average = recode_sev_score_total_tbi_difficulty_concentrating / recode_tbi_difficulty_concentrating_count) %>%
  mutate(recode_tbi_difficulty_remembering_sev_average = recode_sev_score_total_tbi_difficulty_remembering / recode_tbi_difficulty_remembering_count) %>%
  mutate(recode_tbi_blurred_vision_sev_average = recode_sev_score_total_tbi_blurred_vision / recode_tbi_blurred_vision_count) %>%
  mutate(recode_tbi_double_vision_sev_average = recode_sev_score_total_tbi_double_vision / recode_tbi_double_vision_count) %>%
  mutate(recode_tbi_pain_sev_average = recode_sev_score_total_tbi_pain / recode_tbi_pain_count)

######
#total scores for each incidence
test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals %>%
  mutate(recode_qq_covid_1_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_1_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_1_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_1_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_2_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_2_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_2_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_2_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_3_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_3_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_3_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_3_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_4_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_4_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_4_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_4_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_5_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_5_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_5_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_5_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_6_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_6_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_6_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_6_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_7_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_7_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_7_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_7_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_8_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_8_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_8_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_8_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_9_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_9_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_9_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_9_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_10_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_10_sev_score_")))) == ncol(select(., starts_with("recode_qq_covid_10_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_10_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_1_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_1_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_1_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_1_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_2_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_2_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_2_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_2_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_3_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_3_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_3_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_3_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_4_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_4_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_4_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_4_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_5_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_5_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_5_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_5_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_6_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_6_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_6_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_6_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_7_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_7_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_7_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_7_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_8_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_8_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_8_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_8_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_9_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_9_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_9_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_9_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_10_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_10_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_10_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_10_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_11_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_11_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_11_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_11_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_12_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_12_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_12_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_12_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_13_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_13_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_13_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_13_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_14_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_14_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_14_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_14_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_15_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_15_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_15_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_15_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_16_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_16_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_16_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_16_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_17_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_17_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_17_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_17_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_18_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_18_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_18_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_18_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_19_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_19_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_19_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_19_sev_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_20_sev_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_20_sev_score_")))) == ncol(select(., starts_with("recode_qq_tbi_20_sev_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_20_sev_score_")), na.rm=TRUE)))

######
#total scores for each incidence type
test_df_covid_tbi_incidence_totals_types <- test_df_covid_tbi_incidence_totals %>%
  mutate(recode_sev_score_total_covid = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_incidence_total")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_incidence_total"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_sev_score_incidence_total")), na.rm=TRUE))) %>%
  mutate(recode_sev_score_total_tbi = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_incidence_total")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_incidence_total"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_sev_score_incidence_total")), na.rm=TRUE)))

# Total severity score overall
test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals_types %>%
  mutate(recode_overall_sev_score = ifelse(rowSums(is.na(select(., starts_with("recode_sev_score_total_")))) == ncol(select(., starts_with("recode_sev_score_total_"))), NA, rowSums(select(., starts_with("recode_sev_score_total_")), na.rm=TRUE)))
test_df_covid_tbi_incidence_totals$recode_overall_sev_score <- ifelse(test_df_covid_tbi_incidence_totals$recode_overall_sev_score < 0, 0, test_df_covid_tbi_incidence_totals$recode_overall_sev_score)

test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals %>% 
  mutate(covid_symptom_count = rowSums(!is.na(select(., grep("recode_qq_covid_[0-9]+_sev_score_", names(test_df_covid_tbi_incidence_totals)))))) %>%
  mutate(tbi_symptom_count = rowSums(!is.na(select(., grep("recode_qq_tbi_[0-9]+_sev_score_", names(test_df_covid_tbi_incidence_totals)))))) %>%
  mutate(total_symptom_count = covid_symptom_count + tbi_symptom_count)

test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals %>%
  mutate(covid_symptom_sev_average = recode_sev_score_total_covid / covid_symptom_count) %>%
  mutate(tbi_symptom_sev_average = recode_sev_score_total_tbi / tbi_symptom_count) %>%
  mutate(total_symptom_sev_average = recode_overall_sev_score / total_symptom_count)

# Total scores for neurological and psychological symptoms only
#recode_sev_score_total_covid_hypoethesia
covid_neuro_psych_symptom_score_list <- c("recode_sev_score_total_covid_fatigue",
                                          "recode_sev_score_total_covid_insomnia_sleep_problems",
                                          "recode_sev_score_total_covid_mood_swings_irritability",
                                          "recode_sev_score_total_covid_drowsiness",
                                          "recode_sev_score_total_covid_reduced_blurred_vision",
                                          "recode_sev_score_total_covid_photophobia_phonophobia_se",
                                          "recode_sev_score_total_covid_brain_fog",
                                          "recode_sev_score_total_covid_confusion",
                                          "recode_sev_score_total_covid_memory_problems",
                                          "recode_sev_score_total_covid_difficulty_concentrating",
                                          "recode_sev_score_total_covid_delirium",
                                          "recode_sev_score_total_covid_difficulty_finding_words",
                                          "recode_sev_score_total_covid_paresthesia",
                                          "recode_sev_score_total_covid_headache",
                                          "recode_sev_score_total_covid_los",
                                          "recode_sev_score_total_covid_lot",
                                          "recode_sev_score_total_covid_dizziness_lightheadedness",
                                          "recode_sev_score_total_covid_difficulty_balancing",
                                          "recode_sev_score_total_covid_tremors",
                                          "recode_sev_score_total_covid_stroke",
                                          "recode_sev_score_total_covid_seizures",
                                          "recode_sev_score_total_covid_hypoacusis",
                                          "recode_sev_score_total_covid_numbness_hands_feet",
                                          "recode_sev_score_total_covid_hypoethesia")
test_df_covid_tbi_incidence_totals$recode_covid_neuro_psych_sev_score <- rowSums(test_df_covid_tbi_incidence_totals[, covid_neuro_psych_symptom_score_list], na.rm = TRUE)
#recode_sev_score_total_tbi_double_vision
tbi_neuro_psych_symptom_scores_list <- c("recode_sev_score_total_tbi_headache",
                                         "recode_sev_score_total_tbi_balance_problems",
                                         "recode_sev_score_total_tbi_dizziness",
                                         "recode_sev_score_total_tbi_lightheadedness",
                                         "recode_sev_score_total_tbi_fatigue",
                                         "recode_sev_score_total_tbi_trouble_falling_asleep",
                                         "recode_sev_score_total_tbi_sleeping_more",
                                         "recode_sev_score_total_tbi_sleeping_less",
                                         "recode_sev_score_total_tbi_drowsiness",
                                         "recode_sev_score_total_tbi_light_sensitivity",
                                         "recode_sev_score_total_tbi_noise_sensitivity",
                                         "recode_sev_score_total_tbi_irritability",
                                         "recode_sev_score_total_tbi_feeling_frustrated_impatient",
                                         "recode_sev_score_total_tbi_taking_longer_to_think",
                                         "recode_sev_score_total_tbi_restlessness",
                                         "recode_sev_score_total_tbi_sadness",
                                         "recode_sev_score_total_tbi_nervousness_anxiousness",
                                         "recode_sev_score_total_tbi_feeling_more_emotional",
                                         "recode_sev_score_total_tbi_numbness_tingling",
                                         "recode_sev_score_total_tbi_feeling_slowed_down",
                                         "recode_sev_score_total_tbi_in_a_fog",
                                         "recode_sev_score_total_tbi_difficulty_concentrating",
                                         "recode_sev_score_total_tbi_difficulty_remembering",
                                         "recode_sev_score_total_tbi_blurred_vision",
                                         "recode_sev_score_total_tbi_double_vision")

test_df_covid_tbi_incidence_totals$recode_tbi_neuro_psych_sev_score <- rowSums(test_df_covid_tbi_incidence_totals[, tbi_neuro_psych_symptom_scores_list], na.rm = TRUE)
test_df_covid_tbi_incidence_totals$recode_overall_neuro_psych_sev_score <- rowSums(test_df_covid_tbi_incidence_totals[, c("recode_covid_neuro_psych_sev_score", "recode_tbi_neuro_psych_sev_score")], na.rm = TRUE)
test_df_covid_tbi_incidence_totals$recode_overall_neuro_psych_sev_score  <- ifelse(test_df_covid_tbi_incidence_totals$recode_overall_neuro_psych_sev_score  < 0, 0, test_df_covid_tbi_incidence_totals$recode_overall_neuro_psych_sev_score)

test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals %>% 
  mutate(covid_neuro_psych_symptom_count = rowSums(!is.na(select(., all_of(covid_neuro_psych_symptom_score_list))))) %>%
  mutate(tbi_neuro_psych_symptom_count = rowSums(!is.na(select(., all_of(tbi_neuro_psych_symptom_scores_list))))) %>%
  mutate(total_neuro_psych_symptom_count = covid_neuro_psych_symptom_count + tbi_neuro_psych_symptom_count)

test_df_covid_tbi_incidence_totals <- test_df_covid_tbi_incidence_totals %>%
  mutate(covid_neuro_psych_symptom_sev_average = recode_covid_neuro_psych_sev_score / covid_neuro_psych_symptom_count) %>%
  mutate(tbi_neuro_psych_symptom_sev_average = recode_tbi_neuro_psych_sev_score / tbi_neuro_psych_symptom_count) %>%
  mutate(total_neuro_psych_symptom_sev_average = recode_overall_neuro_psych_sev_score / total_neuro_psych_symptom_count)

#####
#age grouping
COAST_Study_Data_Clean_Age_Groups <- test_df_covid_tbi_incidence_totals

COAST_Study_Data_Clean_Age_Groups$Age_Group <- quantcut(COAST_Study_Data_Clean_Age_Groups$age_years, q=2, na.rm=TRUE)
table(COAST_Study_Data_Clean_Age_Groups$Age_Group)

#Need at least 12 young/old to make cell experiments work -  cut out middle 4 middle aged adults 
COAST_Study_Data_Clean_Age_Groups$Age_Group_Long <- ifelse(COAST_Study_Data_Clean_Age_Groups$age_years <= 50.66940, "Young Adults", ifelse(COAST_Study_Data_Clean_Age_Groups$age_years >= 53.56605, "Older Adults","Middle Aged Adults"))
table(COAST_Study_Data_Clean_Age_Groups$Age_Group_Long)

COAST_Study_Data_Clean_Age_Groups$Age_Group_Tert_Long <- ifelse(COAST_Study_Data_Clean_Age_Groups$age_years >=18 & COAST_Study_Data_Clean_Age_Groups$age_years <40, "Young Adults", ifelse(COAST_Study_Data_Clean_Age_Groups$age_years >= 40 & COAST_Study_Data_Clean_Age_Groups$age_years <65, "Middle Aged Adults", ifelse(COAST_Study_Data_Clean_Age_Groups$age_years >= 65 & COAST_Study_Data_Clean_Age_Groups$age_years <= 82, "Older Adults","NA")))
table(COAST_Study_Data_Clean_Age_Groups$Age_Group_Tert_Long)
aggregate(COAST_Study_Data_Clean_Age_Groups$age_years,by=list(COAST_Study_Data_Clean_Age_Groups$Age_Group_Tert_Long),FUN=range,na.rm=TRUE)

#chronic grouping
COAST_Study_Data_Clean_Age_Groups$chronic_covid <- apply(COAST_Study_Data_Clean_Age_Groups %>% select(., starts_with("qq_covid_") & contains("duration")), 1, function(x) any(c("3 months-6 months", "6 months-1 year", ">1 year", "I am still experiencing this symptom") %in% x))
COAST_Study_Data_Clean_Age_Groups$chronic_tbi <- apply(COAST_Study_Data_Clean_Age_Groups %>% select(., starts_with("qq_tbi_") & contains("duration")), 1, function(x) any(c("3 months-6 months", "6 months-1 year", ">1 year", "I am still experiencing this symptom") %in% x))
COAST_Study_Data_Clean_Age_Groups$chronic_acute_overall <- ifelse(COAST_Study_Data_Clean_Age_Groups$chronic_covid == TRUE | COAST_Study_Data_Clean_Age_Groups$chronic_tbi == TRUE, "chronic", "acute")

#####
#frequency
covid__tbi_symptom_freq_columns_to_recode <-
  grep("_freq_",
       names(COAST_Study_Data_Clean_Age_Groups),
       value = TRUE)

#covid__tbi_symptom_freq_columns_to_recode
for (col_name in covid__tbi_symptom_freq_columns_to_recode) {
  COAST_Study_Data_Clean_Age_Groups[paste0("recode_", col_name)] <-
    ifelse(
      COAST_Study_Data_Clean_Age_Groups[col_name] == "Never or almost never had/have symptom",
      1,
      ifelse(
        COAST_Study_Data_Clean_Age_Groups[col_name] == "Sometimes had the symptom",
        2,
        ifelse(
          COAST_Study_Data_Clean_Age_Groups[col_name] == "Often had the symptom",
          3,
          ifelse(
            COAST_Study_Data_Clean_Age_Groups[col_name] == "Frequently had the symptom",
            4,
            ifelse(
              COAST_Study_Data_Clean_Age_Groups[col_name] == "Always had the symptom",
              5,
              ifelse(
                COAST_Study_Data_Clean_Age_Groups[col_name] %in% c(
                  "Never or almost never had/have symptom",
                  "Sometimes had the symptom",
                  "Often had the symptom",
                  "Frequently had the symptom",
                  "Always had the symptom"
                ),
                COAST_Study_Data_Clean_Age_Groups[col_name],
                NA
              )
            )
          )
        )
      )
    )
}

#####
covid_symptom_list_freq <- c("recode_qq_covid_%d_freq_%s_fatigue",
                             "recode_qq_covid_%d_freq_%s_insomnia_sleep_problems",
                             "recode_qq_covid_%d_freq_%s_mood_swings_or_irritability",
                             "recode_qq_covid_%d_freq_%s_drowsiness",
                             "recode_qq_covid_%d_freq_%s_reduced_blurred_vision",
                             "recode_qq_covid_%d_freq_%s_photophobia_phonophobia",
                             "recode_qq_covid_%d_freq_%s_brain_fog",
                             "recode_qq_covid_%d_freq_%s_confusion",
                             "recode_qq_covid_%d_freq_%s_memory_problems",
                             "recode_qq_covid_%d_freq_%s_difficulty_concentrating",
                             "recode_qq_covid_%d_freq_%s_delerium",
                             "recode_qq_covid_%d_freq_%s_difficulty_finding_words",
                             "recode_qq_covid_%d_freq_%s_paresthesia",
                             "recode_qq_covid_%d_freq_%s_headache",
                             "recode_qq_covid_%d_freq_%s_los",
                             "recode_qq_covid_%d_freq_%s_lot",
                             "recode_qq_covid_%d_freq_%s_dizziness_lightheadedness",
                             "recode_qq_covid_%d_freq_%s_difficulty_balancing",
                             "recode_qq_covid_%d_freq_%s_tremors",
                             "recode_qq_covid_%d_freq_%s_stroke",
                             "recode_qq_covid_%d_freq_%s_seizures",
                             "recode_qq_covid_%d_freq_%s_hypoacusis",
                             "recode_qq_covid_%d_freq_%s_numbness_hands_feet",
                             "recode_qq_covid_%d_freq_%s_hypoesthesia",  
                             "recode_qq_covid_%d_freq_%s_nasal_congestion", 
                             "recode_qq_covid_%d_freq_%s_sore_throat", 
                             "recode_qq_covid_%d_freq_%s_runny_nose", 
                             "recode_qq_covid_%d_freq_%s_ear_pain",
                             "recode_qq_covid_%d_freq_%s_cough", 
                             "recode_qq_covid_%d_freq_%s_sputum_production", 
                             "recode_qq_covid_%d_freq_%s_difficulty_breathing_sob", 
                             "recode_qq_covid_%d_freq_%s_hoarse_voice", 
                             "recode_qq_covid_%d_freq_%s_chest_pain_tightness", 
                             "recode_qq_covid_%d_freq_%s_chills", 
                             "recode_qq_covid_%d_freq_%s_swollen_lymph_nodes", 
                             "recode_qq_covid_%d_freq_%s_skipping_meals_appetite_loss", 
                             "recode_qq_covid_%d_freq_%s_sensitivity_heat_cold", 
                             "recode_qq_covid_%d_freq_%s_sweats", 
                             "recode_qq_covid_%d_freq_%s_white_red_purple_swollen_fingers_toes", 
                             "recode_qq_covid_%d_freq_%s_fever_feverish", 
                             "recode_qq_covid_%d_freq_%s_weight_loss", 
                             "recode_qq_covid_%d_freq_%s_tachycardia_arrhythmia_palpitations", 
                             "recode_qq_covid_%d_freq_%s_eye_soreness_discomfort", 
                             "recode_qq_covid_%d_freq_%s_abdominal_pain_stomachache", 
                             "recode_qq_covid_%d_freq_%s_diarrhea", 
                             "recode_qq_covid_%d_freq_%s_nausea_vomiting", 
                             "recode_qq_covid_%d_freq_%s_muscle_weakness", 
                             "recode_qq_covid_%d_freq_%s_muscle_pain_aches", 
                             "recode_qq_covid_%d_freq_%s_bone_and_joint_pain",
                             "recode_qq_covid_%d_freq_%s_neck_back_pain")
tbi_symptom_list_freq <- c("recode_qq_tbi_%d_freq_%s_headache",
                           "recode_qq_tbi_%d_freq_%s_nausea",
                           "recode_qq_tbi_%d_freq_%s_vomiting",
                           "recode_qq_tbi_%d_freq_%s_balance_problems",
                           "recode_qq_tbi_%d_freq_%s_dizziness",
                           "recode_qq_tbi_%d_freq_%s_lightheadedness",
                           "recode_qq_tbi_%d_freq_%s_fatigue",
                           "recode_qq_tbi_%d_freq_%s_trouble_falling_asleep",
                           "recode_qq_tbi_%d_freq_%s_sleeping_more",
                           "recode_qq_tbi_%d_freq_%s_sleeping_less",
                           "recode_qq_tbi_%d_freq_%s_drowsiness",
                           "recode_qq_tbi_%d_freq_%s_light_sensitivity",
                           "recode_qq_tbi_%d_freq_%s_noise_sensitivity",
                           "recode_qq_tbi_%d_freq_%s_irritability",
                           "recode_qq_tbi_%d_freq_%s_feeling_frustrated_impatient",
                           "recode_qq_tbi_%d_freq_%s_taking_longer_to_think",
                           "recode_qq_tbi_%d_freq_%s_restlessness",
                           "recode_qq_tbi_%d_freq_%s_sadness",
                           "recode_qq_tbi_%d_freq_%s_nervousness_anxiousness",
                           "recode_qq_tbi_%d_freq_%s_feeling_more_emotional",
                           "recode_qq_tbi_%d_freq_%s_numbness_tingling",
                           "recode_qq_tbi_%d_freq_%s_feeling_slowed_down",
                           "recode_qq_tbi_%d_freq_%s_in_a_fog",
                           "recode_qq_tbi_%d_freq_%s_difficulty_concentrating",
                           "recode_qq_tbi_%d_freq_%s_difficulty_remembering",
                           "recode_qq_tbi_%d_freq_%s_blurred_vision",
                           "recode_qq_tbi_%d_freq_%s_double_vision", 
                           "recode_qq_tbi_%d_freq_%s_pain")
#####
#rename some columns that are different and are breaking the symptom_incidence_score func
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_1_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_1_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_1_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_2_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_2_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_2_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_3_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_3_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_3_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_4_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_4_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_4_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_5_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_5_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_5_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_6_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_6_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_6_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_7_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_7_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_7_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_8_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_8_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_8_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_9_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_9_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_9_freq_after_mood_swings_irritability <- NULL
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_10_freq_after_mood_swings_or_irritability <- COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_10_freq_after_mood_swings_irritability
COAST_Study_Data_Clean_Age_Groups$recode_qq_covid_10_freq_after_mood_swings_irritability <- NULL

#####
#generate lists of lists
covid_freq_list_of_lists <- generate_list_of_lists(10, covid_symptom_list_freq)
tbi_freq_list_of_lists <- generate_list_of_lists(20, tbi_symptom_list_freq)

#assign new df for freq
COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups

#create single column for each symptom, each incidence
COAST_Study_Data_Clean_Age_Groups_add_freq <- symptom_indicence_score(covid_freq_list_of_lists, COAST_Study_Data_Clean_Age_Groups_add_freq)
COAST_Study_Data_Clean_Age_Groups_add_freq <- symptom_indicence_score(tbi_freq_list_of_lists, COAST_Study_Data_Clean_Age_Groups_add_freq)
#####
#aggregate across symptom
#this works for each symptom, make sure you switch tbi for covid when doing those symptoms and stick to the same naming convention. you will need to change the data frame that plugs nto the next step as well
COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups_add_freq %>%
  mutate(recode_freq_score_total_covid_fatigue = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_fatigue")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_fatigue"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_fatigue")), na.rm=TRUE))) %>%
  mutate(recode_freq_score_total_covid_headache = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_headache")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_headache"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_headache")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_nasal_congestion = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_nasal_congestion")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_nasal_congestion"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_nasal_congestion")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_sore_throat = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sore_throat")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sore_throat"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sore_throat")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_runny_nose = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_runny_nose")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_runny_nose"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_runny_nose")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_ear_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_ear_pain")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_ear_pain"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_ear_pain")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_cough = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_cough")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_cough"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_cough")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_sputum_production = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sputum_production")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sputum_production"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sputum_production")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_difficulty_breathing_sob = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_breathing_sob")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_breathing_sob"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_breathing_sob")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_hoarse_voice = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hoarse_voice")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hoarse_voice"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hoarse_voice")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_chest_pain_tightness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_chest_pain_tightness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_chest_pain_tightness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_chest_pain_tightness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_chills = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_chills")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_chills"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_chills")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_swollen_lymph_nodes = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_swollen_lymph_nodes")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_swollen_lymph_nodes"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_swollen_lymph_nodes")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_skipping_meals_appetite_loss = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_skipping_meals_appetite_loss")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_skipping_meals_appetite_loss"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_skipping_meals_appetite_loss")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_insomnia_sleep_problems = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_insomnia_sleep_problems")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_insomnia_sleep_problems"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_insomnia_sleep_problems")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_sensitivity_heat_cold = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sensitivity_heat_cold")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sensitivity_heat_cold"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sensitivity_heat_cold")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_sweats = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sweats")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sweats"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_sweats")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_white_red_purple_swollen_fingers_toes = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_white_red_purple_swollen_fingers_toes")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_white_red_purple_swollen_fingers_toes"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_white_red_purple_swollen_fingers_toes")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_fever_feverish = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_fever_feverish")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_fever_feverish"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_fever_feverish")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_mood_swings_irritability = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_mood_swings_or_irritability")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_mood_swings_or_irritability"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_mood_swings_or_irritability")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_weight_loss = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_weight_loss")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_weight_loss"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_weight_loss")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_drowsiness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_drowsiness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_drowsiness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_drowsiness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_tachycardia_arrhythmia_palpitations = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_tachycardia_arrhythmia_palpitations")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_tachycardia_arrhythmia_palpitations"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_tachycardia_arrhythmia_palpitations")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_eye_soreness_discomfort = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_eye_soreness_discomfort")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_eye_soreness_discomfort"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_eye_soreness_discomfort")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_reduced_blurred_vision = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_reduced_blurred_vision")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_reduced_blurred_vision"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_reduced_blurred_vision")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_photophobia_phonophobia_se = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_photophobia_phonophobia")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_photophobia_phonophobia"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_photophobia_phonophobia")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_brain_fog = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_brain_fog")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_brain_fog"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_brain_fog")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_confusion = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_confusion")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_confusion"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_confusion")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_memory_problems = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_memory_problems")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_memory_problems"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_memory_problems")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_difficulty_concentrating = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_concentrating")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_concentrating"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_concentrating")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_delirium = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_delirium")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_delirium"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_delirium")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_difficulty_finding_words = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_finding_words")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_finding_words"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_finding_words")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_paresthesia = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_paresthesia")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_paresthesia"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_paresthesia")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_los = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_los")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_los"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_los")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_lot = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_lot")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_lot"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_lot")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_dizziness_lightheadedness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_dizziness_lightheadedness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_dizziness_lightheadedness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_dizziness_lightheadedness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_difficulty_balancing = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_balancing")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_balancing"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_difficulty_balancing")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_tremors = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_tremors")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_tremors"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_tremors")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_stroke = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_stroke")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_stroke"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_stroke")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_seizures = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_seizures")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_seizures"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_seizures")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_hypoacusis = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hypoacusis")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hypoacusis"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hypoacusis")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_numbness_hands_feet = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_numbness_hands_feet")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_numbness_hands_feet"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_numbness_hands_feet")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_hypoethesia = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hypoethesia")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hypoethesia"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_hypoethesia")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_abdominal_pain_stomachache = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_abdominal_pain_stomachache")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_abdominal_pain_stomachache"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_abdominal_pain_stomachache")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_diarrhea = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_diarrhea")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_diarrhea"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_diarrhea")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_nausea_vomiting = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_nausea_vomiting")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_nausea_vomiting"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_nausea_vomiting")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_muscle_weakness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_muscle_weakness")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_muscle_weakness"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_muscle_weakness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_muscle_pain_aches = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_muscle_pain_aches")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_muscle_pain_aches"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_muscle_pain_aches")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_bone_and_joint_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_bone_and_joint_pain")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_bone_and_joint_pain"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_bone_and_joint_pain")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_covid_neck_back_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_neck_back_pain")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_neck_back_pain"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_neck_back_pain")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_headache = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_headache")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_headache"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_headache")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_nausea = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_nausea")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_nausea"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_nausea")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_vomiting = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_vomiting")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_vomiting"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_vomiting")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_balance_problems = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_balance_problems")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_balance_problems"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_balance_problems")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_dizziness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_dizziness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_dizziness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_dizziness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_lightheadedness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_lightheadedness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_lightheadedness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_lightheadedness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_fatigue = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_fatigue")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_fatigue"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_fatigue")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_trouble_falling_asleep = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_trouble_falling_asleep")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_trouble_falling_asleep"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_trouble_falling_asleep")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_sleeping_more = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sleeping_more")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sleeping_more"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sleeping_more")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_sleeping_less = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sleeping_less")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sleeping_less"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sleeping_less")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_drowsiness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_drowsiness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_drowsiness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_drowsiness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_light_sensitivity = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_light_sensitivity")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_light_sensitivity"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_light_sensitivity")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_noise_sensitivity = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_noise_sensitivity")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_noise_sensitivity"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_noise_sensitivity")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_irritability = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_irritability")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_irritability"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_irritability")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_feeling_frustrated_impatient = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_frustrated_impatient")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_frustrated_impatient"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_frustrated_impatient")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_taking_longer_to_think = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_taking_longer_to_think")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_taking_longer_to_think"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_taking_longer_to_think")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_restlessness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_restlessness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_restlessness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_restlessness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_sadness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sadness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sadness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_sadness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_nervousness_anxiousness = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_nervousness_anxiousness")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_nervousness_anxiousness"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_nervousness_anxiousness")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_feeling_more_emotional = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_more_emotional")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_more_emotional"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_more_emotional")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_numbness_tingling = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_numbness_tingling")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_numbness_tingling"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_numbness_tingling")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_feeling_slowed_down = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_slowed_down")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_slowed_down"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_feeling_slowed_down")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_in_a_fog = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_in_a_fog")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_in_a_fog"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_in_a_fog")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_difficulty_concentrating = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_difficulty_concentrating")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_difficulty_concentrating"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_difficulty_concentrating")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_difficulty_remembering = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_difficulty_remembering")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_difficulty_remembering"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_difficulty_remembering")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_blurred_vision = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_blurred_vision")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_blurred_vision"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_blurred_vision")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_double_vision = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_double_vision")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_double_vision"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_double_vision")), na.rm=TRUE))) %>% 
  mutate(recode_freq_score_total_tbi_pain = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_pain")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_pain"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_pain")), na.rm=TRUE))) 

#####
#calculate freq averages
COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups_add_freq %>%
  mutate(recode_covid_fatigue_freq_average = recode_freq_score_total_covid_fatigue / recode_covid_fatigue_count) %>%
  mutate(recode_covid_headache_freq_average = recode_freq_score_total_covid_headache / recode_covid_headache_count) %>%
  mutate(recode_covid_nasal_congestion_freq_average = recode_freq_score_total_covid_nasal_congestion / recode_covid_nasal_congestion_count) %>%
  mutate(recode_covid_sore_throat_freq_average = recode_freq_score_total_covid_sore_throat / recode_covid_sore_throat_count) %>%
  mutate(recode_covid_runny_nose_freq_average = recode_freq_score_total_covid_runny_nose / recode_covid_runny_nose_count) %>%
  mutate(recode_covid_ear_pain_freq_average = recode_freq_score_total_covid_ear_pain / recode_covid_ear_pain_count) %>%
  mutate(recode_covid_cough_freq_average = recode_freq_score_total_covid_cough / recode_covid_cough_count) %>%
  mutate(recode_covid_sputum_production_freq_average = recode_freq_score_total_covid_sputum_production / recode_covid_sputum_production_count) %>%
  mutate(recode_covid_difficulty_breathing_sob_freq_average = recode_freq_score_total_covid_difficulty_breathing_sob / recode_covid_difficulty_breathing_sob_count) %>%
  mutate(recode_covid_hoarse_voice_freq_average = recode_freq_score_total_covid_hoarse_voice / recode_covid_hoarse_voice_count) %>%
  mutate(recode_covid_chest_pain_tightness_freq_average = recode_freq_score_total_covid_chest_pain_tightness / recode_covid_chest_pain_tightness_count) %>%
  mutate(recode_covid_chills_freq_average = recode_freq_score_total_covid_chills / recode_covid_chills_count) %>%
  mutate(recode_covid_swollen_lymph_nodes_freq_average = recode_freq_score_total_covid_swollen_lymph_nodes / recode_covid_swollen_lymph_nodes_count) %>%
  mutate(recode_covid_skipping_meals_appetite_loss_freq_average = recode_freq_score_total_covid_skipping_meals_appetite_loss / recode_covid_skipping_meals_appetite_loss_count) %>%
  mutate(recode_covid_insomnia_sleep_problems_freq_average = recode_freq_score_total_covid_insomnia_sleep_problems / recode_covid_insomnia_sleep_problems_count) %>%
  mutate(recode_covid_sensitivity_heat_cold_freq_average = recode_freq_score_total_covid_sensitivity_heat_cold / recode_covid_sensitivity_heat_cold_count) %>%
  mutate(recode_covid_sweats_freq_average = recode_freq_score_total_covid_sweats / recode_covid_sweats_count) %>%
  mutate(recode_covid_white_red_purple_swollen_fingers_toes_freq_average = recode_freq_score_total_covid_white_red_purple_swollen_fingers_toes / recode_covid_white_red_purple_swollen_fingers_toes_count) %>%
  mutate(recode_covid_fever_feverish_freq_average = recode_freq_score_total_covid_fever_feverish / recode_covid_fever_feverish_count) %>%
  mutate(recode_covid_mood_swings_irritability_freq_average = recode_freq_score_total_covid_mood_swings_irritability / recode_covid_mood_swings_irritability_count) %>%
  mutate(recode_covid_weight_loss_freq_average = recode_freq_score_total_covid_weight_loss / recode_covid_weight_loss_count) %>%
  mutate(recode_covid_drowsiness_freq_average = recode_freq_score_total_covid_drowsiness / recode_covid_drowsiness_count) %>%
  mutate(recode_covid_tachycardia_arrhythmia_palpitations_freq_average = recode_freq_score_total_covid_tachycardia_arrhythmia_palpitations / recode_covid_tachycardia_arrhythmia_palpitations_count) %>%
  mutate(recode_covid_eye_soreness_discomfort_freq_average = recode_freq_score_total_covid_eye_soreness_discomfort / recode_covid_eye_soreness_discomfort_count) %>%
  mutate(recode_covid_reduced_blurred_vision_freq_average = recode_freq_score_total_covid_reduced_blurred_vision / recode_covid_reduced_blurred_vision_count) %>%
  mutate(recode_covid_photophobia_phonophobia_se_freq_average = recode_freq_score_total_covid_photophobia_phonophobia_se / recode_covid_photophobia_phonophobia_se_count) %>%
  mutate(recode_covid_brain_fog_freq_average = recode_freq_score_total_covid_brain_fog / recode_covid_brain_fog_count) %>%
  mutate(recode_covid_confusion_freq_average = recode_freq_score_total_covid_confusion / recode_covid_confusion_count) %>%
  mutate(recode_covid_memory_problems_freq_average = recode_freq_score_total_covid_memory_problems / recode_covid_memory_problems_count) %>%
  mutate(recode_covid_difficulty_concentrating_freq_average = recode_freq_score_total_covid_difficulty_concentrating / recode_covid_difficulty_concentrating_count) %>%
  mutate(recode_covid_delirium_freq_average = recode_freq_score_total_covid_delirium / recode_covid_delirium_count) %>%
  mutate(recode_covid_difficulty_finding_words_freq_average = recode_freq_score_total_covid_difficulty_finding_words / recode_covid_difficulty_finding_words_count) %>%
  mutate(recode_covid_paresthesia_freq_average = recode_freq_score_total_covid_paresthesia / recode_covid_paresthesia_count) %>%
  mutate(recode_covid_los_freq_average = recode_freq_score_total_covid_los / recode_covid_los_count) %>%
  mutate(recode_covid_lot_freq_average = recode_freq_score_total_covid_lot / recode_covid_lot_count) %>%
  mutate(recode_covid_dizziness_lightheadedness_freq_average = recode_freq_score_total_covid_dizziness_lightheadedness / recode_covid_dizziness_lightheadedness_count) %>%
  mutate(recode_covid_difficulty_balancing_freq_average = recode_freq_score_total_covid_difficulty_balancing / recode_covid_difficulty_balancing_count) %>%
  mutate(recode_covid_tremors_freq_average = recode_freq_score_total_covid_tremors / recode_covid_tremors_count) %>%
  mutate(recode_covid_stroke_freq_average = recode_freq_score_total_covid_stroke / recode_covid_stroke_count) %>%
  mutate(recode_covid_seizures_freq_average = recode_freq_score_total_covid_seizures / recode_covid_seizures_count) %>%
  mutate(recode_covid_hypoacusis_freq_average = recode_freq_score_total_covid_hypoacusis / recode_covid_hypoacusis_count) %>%
  mutate(recode_covid_numbness_hands_feet_freq_average = recode_freq_score_total_covid_numbness_hands_feet / recode_covid_numbness_hands_feet_count) %>%
  mutate(recode_covid_hypoethesia_freq_average = recode_freq_score_total_covid_hypoethesia / recode_covid_hypoethesia_count) %>%
  mutate(recode_covid_abdominal_pain_stomachache_freq_average = recode_freq_score_total_covid_abdominal_pain_stomachache / recode_covid_abdominal_pain_stomachache_count) %>%
  mutate(recode_covid_diarrhea_freq_average = recode_freq_score_total_covid_diarrhea / recode_covid_diarrhea_count) %>%
  mutate(recode_covid_nausea_vomiting_freq_average = recode_freq_score_total_covid_nausea_vomiting / recode_covid_nausea_vomiting_count) %>%
  mutate(recode_covid_muscle_weakness_freq_average = recode_freq_score_total_covid_muscle_weakness / recode_covid_muscle_weakness_count) %>%
  mutate(recode_covid_muscle_pain_aches_freq_average = recode_freq_score_total_covid_muscle_pain_aches / recode_covid_muscle_pain_aches_count) %>%
  mutate(recode_covid_bone_and_joint_pain_freq_average = recode_freq_score_total_covid_bone_and_joint_pain / recode_covid_bone_and_joint_pain_count) %>%
  mutate(recode_covid_neck_back_pain_freq_average = recode_freq_score_total_covid_neck_back_pain / recode_covid_neck_back_pain_count) %>%
  mutate(recode_tbi_headache_freq_average = recode_freq_score_total_tbi_headache / recode_tbi_headache_count) %>%
  mutate(recode_tbi_nausea_freq_average = recode_freq_score_total_tbi_nausea / recode_tbi_nausea_count) %>%
  mutate(recode_tbi_vomiting_freq_average = recode_freq_score_total_tbi_vomiting / recode_tbi_vomiting_count) %>%
  mutate(recode_tbi_balance_problems_freq_average = recode_freq_score_total_tbi_balance_problems / recode_tbi_balance_problems_count) %>%
  mutate(recode_tbi_dizziness_freq_average = recode_freq_score_total_tbi_dizziness / recode_tbi_dizziness_count) %>%
  mutate(recode_tbi_lightheadedness_freq_average = recode_freq_score_total_tbi_lightheadedness / recode_tbi_lightheadedness_count) %>%
  mutate(recode_tbi_fatigue_freq_average = recode_freq_score_total_tbi_fatigue / recode_tbi_fatigue_count) %>%
  mutate(recode_tbi_trouble_falling_asleep_freq_average = recode_freq_score_total_tbi_trouble_falling_asleep / recode_tbi_trouble_falling_asleep_count) %>%
  mutate(recode_tbi_sleeping_more_freq_average = recode_freq_score_total_tbi_sleeping_more / recode_tbi_sleeping_more_count) %>%
  mutate(recode_tbi_sleeping_less_freq_average = recode_freq_score_total_tbi_sleeping_less / recode_tbi_sleeping_less_count) %>%
  mutate(recode_tbi_drowsiness_freq_average = recode_freq_score_total_tbi_drowsiness / recode_tbi_drowsiness_count) %>%
  mutate(recode_tbi_light_sensitivity_freq_average = recode_freq_score_total_tbi_light_sensitivity / recode_tbi_light_sensitivity_count) %>%
  mutate(recode_tbi_noise_sensitivity_freq_average = recode_freq_score_total_tbi_noise_sensitivity / recode_tbi_noise_sensitivity_count) %>%
  mutate(recode_tbi_irritability_freq_average = recode_freq_score_total_tbi_irritability / recode_tbi_irritability_count) %>%
  mutate(recode_tbi_feeling_frustrated_impatient_freq_average = recode_freq_score_total_tbi_feeling_frustrated_impatient / recode_tbi_feeling_frustrated_impatient_count) %>%
  mutate(recode_tbi_taking_longer_to_think_freq_average = recode_freq_score_total_tbi_taking_longer_to_think / recode_tbi_taking_longer_to_think_count) %>%
  mutate(recode_tbi_restlessness_freq_average = recode_freq_score_total_tbi_restlessness / recode_tbi_restlessness_count) %>%
  mutate(recode_tbi_sadness_freq_average = recode_freq_score_total_tbi_sadness / recode_tbi_sadness_count) %>%
  mutate(recode_tbi_nervousness_anxiousness_freq_average = recode_freq_score_total_tbi_nervousness_anxiousness / recode_tbi_nervousness_anxiousness_count) %>%
  mutate(recode_tbi_feeling_more_emotional_freq_average = recode_freq_score_total_tbi_feeling_more_emotional / recode_tbi_feeling_more_emotional_count) %>%
  mutate(recode_tbi_numbness_tingling_freq_average = recode_freq_score_total_tbi_numbness_tingling / recode_tbi_numbness_tingling_count) %>%
  mutate(recode_tbi_feeling_slowed_down_freq_average = recode_freq_score_total_tbi_feeling_slowed_down / recode_tbi_feeling_slowed_down_count) %>%
  mutate(recode_tbi_in_a_fog_freq_average = recode_freq_score_total_tbi_in_a_fog / recode_tbi_in_a_fog_count) %>%
  mutate(recode_tbi_difficulty_concentrating_freq_average = recode_freq_score_total_tbi_difficulty_concentrating / recode_tbi_difficulty_concentrating_count) %>%
  mutate(recode_tbi_difficulty_remembering_freq_average = recode_freq_score_total_tbi_difficulty_remembering / recode_tbi_difficulty_remembering_count) %>%
  mutate(recode_tbi_blurred_vision_freq_average = recode_freq_score_total_tbi_blurred_vision / recode_tbi_blurred_vision_count) %>%
  mutate(recode_tbi_double_vision_freq_average = recode_freq_score_total_tbi_double_vision / recode_tbi_double_vision_count) %>%
  mutate(recode_tbi_pain_freq_average = recode_freq_score_total_tbi_pain / recode_tbi_pain_count)
#total scores for each incidence
COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups_add_freq %>%
  mutate(recode_qq_covid_1_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_1_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_1_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_1_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_2_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_2_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_2_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_2_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_3_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_3_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_3_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_3_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_4_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_4_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_4_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_4_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_5_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_5_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_5_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_5_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_6_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_6_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_6_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_6_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_7_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_7_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_7_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_7_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_8_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_8_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_8_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_8_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_9_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_9_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_9_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_9_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_covid_10_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_10_freq_score_")))) == ncol(select(., starts_with("recode_qq_covid_10_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_covid_10_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_1_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_1_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_1_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_1_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_2_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_2_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_2_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_2_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_3_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_3_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_3_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_3_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_4_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_4_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_4_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_4_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_5_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_5_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_5_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_5_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_6_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_6_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_6_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_6_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_7_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_7_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_7_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_7_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_8_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_8_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_8_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_8_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_9_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_9_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_9_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_9_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_10_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_10_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_10_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_10_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_11_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_11_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_11_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_11_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_12_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_12_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_12_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_12_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_13_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_13_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_13_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_13_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_14_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_14_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_14_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_14_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_15_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_15_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_15_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_15_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_16_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_16_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_16_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_16_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_17_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_17_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_17_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_17_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_18_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_18_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_18_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_18_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_19_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_19_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_19_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_19_freq_score_")), na.rm=TRUE))) %>%
  mutate(recode_qq_tbi_20_freq_score_incidence_total = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_20_freq_score_")))) == ncol(select(., starts_with("recode_qq_tbi_20_freq_score_"))), NA, rowSums(select(., starts_with("recode_qq_tbi_20_freq_score_")), na.rm=TRUE)))

#total scores for each incidence type
COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups_add_freq %>%
  mutate(recode_freq_score_total_covid = ifelse(rowSums(is.na(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_incidence_total")))) == ncol(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_incidence_total"))), NA, rowSums(select(., starts_with("recode_qq_covid_") & ends_with("_freq_score_incidence_total")), na.rm=TRUE))) %>%
  mutate(recode_freq_score_total_tbi = ifelse(rowSums(is.na(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_incidence_total")))) == ncol(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_incidence_total"))), NA, rowSums(select(., starts_with("recode_qq_tbi_") & ends_with("_freq_score_incidence_total")), na.rm=TRUE)))
# COAST_Study_Data_Clean_Age_Groups_add_freq$recode_freq_score_total_covid
# COAST_Study_Data_Clean_Age_Groups_add_freq$recode_freq_score_total_tbi
#total score overall
COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups_add_freq %>%
  mutate(recode_overall_freq_score = ifelse(rowSums(is.na(select(., starts_with("recode_freq_score_total_")))) == ncol(select(., starts_with("recode_freq_score_total_"))), NA, rowSums(select(., starts_with("recode_freq_score_total_")), na.rm=TRUE)))
COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_freq_score <- ifelse(COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_freq_score < 0, 0, COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_freq_score)
# COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_freq_score
# average scores
COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups_add_freq %>%
  mutate(covid_symptom_freq_average = recode_freq_score_total_covid / covid_symptom_count) %>%
  mutate(tbi_symptom_freq_average = recode_freq_score_total_tbi / tbi_symptom_count) %>%
  mutate(total_symptom_freq_average = recode_overall_freq_score / total_symptom_count)

# total scores for neurological and psychological only
#recode_freq_score_total_covid_hypoethesia
covid_neuro_psych_symptom_score_list <- c("recode_freq_score_total_covid_fatigue",
                                          "recode_freq_score_total_covid_insomnia_sleep_problems",
                                          "recode_freq_score_total_covid_mood_swings_irritability",
                                          "recode_freq_score_total_covid_drowsiness",
                                          "recode_freq_score_total_covid_reduced_blurred_vision",
                                          "recode_freq_score_total_covid_photophobia_phonophobia_se",
                                          "recode_freq_score_total_covid_brain_fog",
                                          "recode_freq_score_total_covid_confusion",
                                          "recode_freq_score_total_covid_memory_problems",
                                          "recode_freq_score_total_covid_difficulty_concentrating",
                                          "recode_freq_score_total_covid_delirium",
                                          "recode_freq_score_total_covid_difficulty_finding_words",
                                          "recode_freq_score_total_covid_paresthesia",
                                          "recode_freq_score_total_covid_headache",
                                          "recode_freq_score_total_covid_los",
                                          "recode_freq_score_total_covid_lot",
                                          "recode_freq_score_total_covid_dizziness_lightheadedness",
                                          "recode_freq_score_total_covid_difficulty_balancing",
                                          "recode_freq_score_total_covid_tremors",
                                          "recode_freq_score_total_covid_stroke",
                                          "recode_freq_score_total_covid_seizures",
                                          "recode_freq_score_total_covid_hypoacusis",
                                          "recode_freq_score_total_covid_numbness_hands_feet",
                                          "recode_freq_score_total_covid_hypoethesia")
COAST_Study_Data_Clean_Age_Groups_add_freq$recode_covid_neuro_psych_freq_score <- rowSums(COAST_Study_Data_Clean_Age_Groups_add_freq[, covid_neuro_psych_symptom_score_list], na.rm = TRUE)
#recode_freq_score_total_tbi_double_vision
tbi_neuro_psych_symptom_scores_list <- c("recode_freq_score_total_tbi_headache",
                                         "recode_freq_score_total_tbi_balance_problems",
                                         "recode_freq_score_total_tbi_dizziness",
                                         "recode_freq_score_total_tbi_lightheadedness",
                                         "recode_freq_score_total_tbi_fatigue",
                                         "recode_freq_score_total_tbi_trouble_falling_asleep",
                                         "recode_freq_score_total_tbi_sleeping_more",
                                         "recode_freq_score_total_tbi_sleeping_less",
                                         "recode_freq_score_total_tbi_drowsiness",
                                         "recode_freq_score_total_tbi_light_sensitivity",
                                         "recode_freq_score_total_tbi_noise_sensitivity",
                                         "recode_freq_score_total_tbi_irritability",
                                         "recode_freq_score_total_tbi_feeling_frustrated_impatient",
                                         "recode_freq_score_total_tbi_taking_longer_to_think",
                                         "recode_freq_score_total_tbi_restlessness",
                                         "recode_freq_score_total_tbi_sadness",
                                         "recode_freq_score_total_tbi_nervousness_anxiousness",
                                         "recode_freq_score_total_tbi_feeling_more_emotional",
                                         "recode_freq_score_total_tbi_numbness_tingling",
                                         "recode_freq_score_total_tbi_feeling_slowed_down",
                                         "recode_freq_score_total_tbi_in_a_fog",
                                         "recode_freq_score_total_tbi_difficulty_concentrating",
                                         "recode_freq_score_total_tbi_difficulty_remembering",
                                         "recode_freq_score_total_tbi_blurred_vision",
                                         "recode_freq_score_total_tbi_double_vision")
COAST_Study_Data_Clean_Age_Groups_add_freq$recode_tbi_neuro_psych_freq_score <- rowSums(COAST_Study_Data_Clean_Age_Groups_add_freq[, tbi_neuro_psych_symptom_scores_list], na.rm = TRUE)
COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_neuro_psych_freq_score <- rowSums(COAST_Study_Data_Clean_Age_Groups_add_freq[, c("recode_covid_neuro_psych_freq_score", "recode_tbi_neuro_psych_freq_score")], na.rm = TRUE)
COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_neuro_psych_freq_score <- ifelse(COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_neuro_psych_freq_score < 0, 0, COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_neuro_psych_freq_score)

COAST_Study_Data_Clean_Age_Groups_add_freq <- COAST_Study_Data_Clean_Age_Groups_add_freq %>%
  mutate(covid_neuro_psych_symptom_freq_average = recode_covid_neuro_psych_freq_score / covid_neuro_psych_symptom_count) %>%
  mutate(tbi_neuro_psych_symptom_freq_average = recode_tbi_neuro_psych_freq_score / tbi_neuro_psych_symptom_count) %>%
  mutate(total_neuro_psych_symptom_freq_average = recode_overall_neuro_psych_freq_score / total_neuro_psych_symptom_count)

# Combined severity + frequency score for ranking
COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_sev_plus_freq_score <- COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_freq_score + COAST_Study_Data_Clean_Age_Groups_add_freq$recode_overall_sev_score

#####
test_stats_df <-
  COAST_Study_Data_Clean_Age_Groups_add_freq %>% select(
    "participant_id",
    "Age_Group_Long",
    "age_years",
    "qq_group",
    "qq_biological_sex",
    "chronic_acute_overall",
    "recode_sev_score_total_covid",
    "recode_sev_score_total_tbi",
    "recode_freq_score_total_covid",
    "recode_freq_score_total_tbi",
    "recode_overall_sev_score",
    "recode_overall_freq_score",
    "recode_overall_sev_plus_freq_score",
    "recode_tbi_neuro_psych_sev_score",
    "recode_tbi_neuro_psych_freq_score",
    "recode_covid_neuro_psych_sev_score",
    "recode_covid_neuro_psych_freq_score",
    "recode_overall_neuro_psych_sev_score",
    "recode_overall_neuro_psych_freq_score"
  )

top_5_chronic_young <- test_stats_df %>% filter(Age_Group_Long == "Young Adults" & chronic_acute_overall == "chronic") %>% arrange(-recode_overall_neuro_psych_sev_score) %>% slice(1:5) %>% select("participant_id", "Age_Group_Long", "recode_overall_neuro_psych_sev_score", "age_years")
bot_3_chronic_young_test_stats_df_sev <- test_stats_df %>% filter(Age_Group_Long == "Young Adults" & chronic_acute_overall == "chronic") %>% arrange(recode_overall_sev_score) %>% slice(1:3) %>% select("participant_id", "Age_Group_Long", "recode_overall_sev_score", "age_years")
top_3_acute_young_test_stats_df_sev <- test_stats_df %>% filter(Age_Group_Long == "Young Adults" & chronic_acute_overall == "acute") %>% arrange(-recode_overall_sev_score) %>% slice(1:3) %>% select("participant_id", "Age_Group_Long", "recode_overall_sev_score", "age_years")
bot_3_acute_young_test_stats_df_sev <- test_stats_df %>% filter(Age_Group_Long == "Young Adults" & chronic_acute_overall == "acute") %>% arrange(recode_overall_sev_score) %>% slice(1:3) %>% select("participant_id", "Age_Group_Long", "recode_overall_sev_score", "age_years")

top_3_chronic_old_test_stats_df_sev <- test_stats_df %>% filter(Age_Group_Long == "Older Adults" & chronic_acute_overall == "chronic") %>% arrange(-recode_overall_sev_score) %>% slice(1:3) %>% select("participant_id", "Age_Group_Long", "recode_overall_sev_score", "age_years")
bot_3_chronic_old_test_stats_df_sev <- test_stats_df %>% filter(Age_Group_Long == "Older Adults" & chronic_acute_overall == "chronic") %>% arrange(recode_overall_sev_score) %>% slice(1:3) %>% select("participant_id", "Age_Group_Long", "recode_overall_sev_score", "age_years")
top_3_acute_old_test_stats_df_sev <- test_stats_df %>% filter(Age_Group_Long == "Older Adults" & chronic_acute_overall == "acute") %>% arrange(-recode_overall_sev_score) %>% slice(1:3) %>% select("participant_id", "Age_Group_Long", "recode_overall_sev_score", "age_years")
bot_3_acute_old_test_stats_df_sev <- test_stats_df %>% filter(Age_Group_Long == "Older Adults" & chronic_acute_overall == "acute") %>% arrange(recode_overall_sev_score) %>% slice(1:3) %>% select("participant_id", "Age_Group_Long", "recode_overall_sev_score", "age_years")

table(COAST_Study_Data_Clean_only_covid_pos$qq_group)
controls <- COAST_Study_Data_Clean_only_covid_pos %>% filter(qq_group == "COVID-19 (-) mTBI (-)") %>% select("participant_id", "age_years")

# top 12 highest neuro/psych severity score each group (4 main study groups)
# ID, group, age, gender, score
# top_12_double_pos <- test_stats_df %>% filter(qq_group == "COVID-19 (+) mTBI (+)") %>% arrange(-recode_overall_neuro_psych_sev_score) %>% slice(1:12) %>% select("participant_id", "qq_group", "age_years", "qq_biological_sex", "recode_overall_neuro_psych_sev_score")
# top_12_tbi_only <- test_stats_df %>% filter(qq_group == "COVID-19 (-) mTBI (+)") %>% arrange(-recode_overall_neuro_psych_sev_score) %>% slice(1:12) %>% select("participant_id", "qq_group", "age_years", "qq_biological_sex", "recode_overall_neuro_psych_sev_score")
# top_12_covid_only <- test_stats_df %>% filter(qq_group == "COVID-19 (+) mTBI (-)") %>% arrange(-recode_overall_neuro_psych_sev_score) %>% slice(1:12) %>% select("participant_id", "qq_group", "age_years", "qq_biological_sex", "recode_overall_neuro_psych_sev_score")
# all_controls <- test_stats_df %>% filter(qq_group == "COVID-19 (-) mTBI (-)") %>% arrange(-recode_overall_neuro_psych_sev_score) %>% select("participant_id", "qq_group", "age_years", "qq_biological_sex", "recode_overall_neuro_psych_sev_score")
# (optional: write top-12/controls CSVs to a path of your choice)


#####
#dates
#time between last illness or injury and the study
#date formatting
COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_freq
#covid dates
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_1_pos_test_year <- COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_1_pos_test_year + 1899
COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_covid_1_pos_test_date", qq_covid_1_pos_test_month:qq_covid_1_pos_test_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_1_pos_test_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_1_pos_test_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_1_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_1_pos_test_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')


COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_2_pos_test_year <- COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_2_pos_test_year + 1899
COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_covid_2_pos_test_date", qq_covid_2_pos_test_month:qq_covid_2_pos_test_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_2_pos_test_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_2_pos_test_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_2_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_covid_2_pos_test_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

#tbi dates
TBI_month_columns <- grep("qq_tbi_[0-9]+_month", names(COAST_Study_Data_Clean_Age_Groups_add_dates), value = TRUE)
# TBI_month_columns

for (col_name in TBI_month_columns) {
  COAST_Study_Data_Clean_Age_Groups_add_dates[paste0(col_name)] <- ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "January", "01-15", 
                                                                          ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "February", "02-15", 
                                                                                 ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "March", "03-15", 
                                                                                        ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "April", "04-15", 
                                                                                               ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "May", "05-15", 
                                                                                                      ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "June", "06-15", 
                                                                                                             ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "July", "07-15", 
                                                                                                                    ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "August", "08-15", 
                                                                                                                           ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "September", "09-15", 
                                                                                                                                  ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "October", "10-15", 
                                                                                                                                         ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "November", "11-15", 
                                                                                                                                                ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] == "December", "12-15", 
                                                                                                                                                       ifelse(COAST_Study_Data_Clean_Age_Groups_add_dates[col_name] %in% c("January", "February", "March", "April","May", "June", "July", "August", "September", "October", "November", "December"), 
                                                                                                                                                              COAST_Study_Data_Clean[col_name],NA)))))))))))))
}
# COAST_Study_Data_Clean_Age_Groups_add_dates["qq_tbi_1_month"]

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_1_symptom_onset_date", qq_tbi_1_month, qq_tbi_1_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_1_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_1_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_1_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_1_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_2_symptom_onset_date", qq_tbi_2_month, qq_tbi_2_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_2_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_2_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_2_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_2_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_3_symptom_onset_date", qq_tbi_3_month, qq_tbi_3_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_3_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_3_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_3_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_3_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_4_symptom_onset_date", qq_tbi_4_month, qq_tbi_4_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_4_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_4_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_4_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_4_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_5_symptom_onset_date", qq_tbi_5_month, qq_tbi_5_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_5_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_5_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_5_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_5_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_6_symptom_onset_date", qq_tbi_6_month, qq_tbi_6_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_6_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_6_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_6_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_6_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_7_symptom_onset_date", qq_tbi_7_month, qq_tbi_7_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_7_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_7_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_7_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_7_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_8_symptom_onset_date", qq_tbi_8_month, qq_tbi_8_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_8_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_8_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_8_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_8_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_9_symptom_onset_date", qq_tbi_9_month, qq_tbi_9_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_9_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_9_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_9_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_9_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_10_symptom_onset_date", qq_tbi_10_month, qq_tbi_10_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_10_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_10_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_10_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_10_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_11_symptom_onset_date", qq_tbi_11_month, qq_tbi_11_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_11_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_11_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_11_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_11_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_12_symptom_onset_date", qq_tbi_12_month, qq_tbi_12_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_12_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_12_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_12_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_12_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_13_symptom_onset_date", qq_tbi_13_month, qq_tbi_13_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_13_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_13_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_13_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_13_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_14_symptom_onset_date", qq_tbi_14_month, qq_tbi_14_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_14_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_14_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_14_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_14_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_15_symptom_onset_date", qq_tbi_15_month, qq_tbi_15_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_15_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_15_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_15_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_15_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_16_symptom_onset_date", qq_tbi_16_month, qq_tbi_16_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_16_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_16_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_16_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_16_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_17_symptom_onset_date", qq_tbi_17_month, qq_tbi_17_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_17_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_17_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_17_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_17_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_18_symptom_onset_date", qq_tbi_18_month, qq_tbi_18_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_18_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_18_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_18_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_18_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_19_symptom_onset_date", qq_tbi_19_month, qq_tbi_19_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_19_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_19_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_19_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_19_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% unite("qq_tbi_20_symptom_onset_date", qq_tbi_20_month, qq_tbi_20_year, sep= "-", remove = FALSE)
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_20_symptom_onset_date <- as.Date(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_20_symptom_onset_date,format = "%m-%d-%Y")
COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_20_years_ago <- interval(as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_20_symptom_onset_date), as.POSIXct(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_start_date)) %>% as.numeric('years')

#average time since covid
covid_years_ago_cols <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% select(grep("qq_covid_[0-9]+_years_ago", names(COAST_Study_Data_Clean_Age_Groups_add_dates)))
COAST_Study_Data_Clean_Age_Groups_add_dates$average_years_since_covid <- rowMeans(covid_years_ago_cols, na.rm=T)
#average time since tbi
tbi_years_ago_cols <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% select(grep("qq_tbi_[0-9]+_years_ago", names(COAST_Study_Data_Clean_Age_Groups_add_dates)))
COAST_Study_Data_Clean_Age_Groups_add_dates$average_years_since_tbi <- rowMeans(tbi_years_ago_cols, na.rm=T)

#real number of covid incidences
covid_test_results_cols <- names(COAST_Study_Data_Clean_Age_Groups_add_dates)[grep("qq_covid_[^_]+_test_results", names(COAST_Study_Data_Clean_Age_Groups_add_dates))]
yes_count <- rowSums(COAST_Study_Data_Clean_Age_Groups_add_dates[, covid_test_results_cols] == "Yes", na.rm = TRUE)
COAST_Study_Data_Clean_Age_Groups_add_dates$covid_pos_test_num <- yes_count

##
#File merge 

join_and_select <- function(df_A, df_B, columns_to_select) {
  
  # Ensure 'participant_id' is in both dataframes
  if (!("participant_id" %in% names(df_A)) || !("participant_id" %in% names(df_B))) {
    stop("Both dataframes must have a 'participant_id' column")
  }
  
  # Ensure all columns_to_select are in df_A
  if (!all(columns_to_select %in% names(df_A))) {
    missing_cols <- setdiff(columns_to_select, names(df_A))
    stop(paste("The following columns are not in df_A:", paste(missing_cols, collapse = ", ")))
  }
  
  # Select specified columns from df_A and join with df_B
  result <- df_A %>%
    select(participant_id, all_of(columns_to_select)) %>%
    right_join(df_B, by = "participant_id")
  
  return(result)
}
# total incidence nums for all incidence types
COAST_Study_Data_Clean_Age_Groups_add_dates <- COAST_Study_Data_Clean_Age_Groups_add_dates %>%
  mutate(qq_tbi_num = replace_na(qq_tbi_num, 0))
COAST_Study_Data_Clean_Age_Groups_add_dates$total_incidences <- COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_num + COAST_Study_Data_Clean_Age_Groups_add_dates$covid_pos_test_num

# calculate average number of symptoms reported per incidence
COAST_Study_Data_Clean_Age_Groups_add_dates$covid_total_symptoms_reported <- rowSums(COAST_Study_Data_Clean_Age_Groups_add_dates[, c(
  "recode_covid_fatigue_count",
  "recode_covid_headache_count",
  "recode_covid_insomnia_sleep_problems_count",
  "recode_covid_mood_swings_irritability_count",
  "recode_covid_drowsiness_count",
  "recode_covid_reduced_blurred_vision_count",
  "recode_covid_photophobia_phonophobia_se_count",
  "recode_covid_brain_fog_count",
  "recode_covid_confusion_count",
  "recode_covid_memory_problems_count",
  "recode_covid_difficulty_concentrating_count",
  "recode_covid_delirium_count",
  "recode_covid_difficulty_finding_words_count",
  "recode_covid_paresthesia_count",
  "recode_covid_los_count",
  "recode_covid_lot_count",
  "recode_covid_dizziness_lightheadedness_count",
  "recode_covid_difficulty_balancing_count",
  "recode_covid_tremors_count",
  "recode_covid_stroke_count",
  "recode_covid_seizures_count",
  "recode_covid_hypoacusis_count",
  "recode_covid_numbness_hands_feet_count",
  "recode_covid_hypoethesia_count"
)], na.rm = TRUE)

COAST_Study_Data_Clean_Age_Groups_add_dates$tbi_total_symptoms_reported <- rowSums(COAST_Study_Data_Clean_Age_Groups_add_dates[, c(
  "recode_tbi_headache_count",
  "recode_tbi_balance_problems_count",
  "recode_tbi_dizziness_count",
  "recode_tbi_lightheadedness_count",
  "recode_tbi_fatigue_count",
  "recode_tbi_trouble_falling_asleep_count",
  "recode_tbi_sleeping_more_count",
  "recode_tbi_sleeping_less_count",
  "recode_tbi_drowsiness_count",
  "recode_tbi_light_sensitivity_count",
  "recode_tbi_noise_sensitivity_count",
  "recode_tbi_irritability_count",
  "recode_tbi_feeling_frustrated_impatient_count",
  "recode_tbi_taking_longer_to_think_count",
  "recode_tbi_restlessness_count",
  "recode_tbi_sadness_count",
  "recode_tbi_nervousness_anxiousness_count",
  "recode_tbi_feeling_more_emotional_count",
  "recode_tbi_numbness_tingling_count",
  "recode_tbi_feeling_slowed_down_count",
  "recode_tbi_in_a_fog_count",
  "recode_tbi_difficulty_concentrating_count",
  "recode_tbi_difficulty_remembering_count",
  "recode_tbi_blurred_vision_count",
  "recode_tbi_double_vision_count"
)], na.rm = TRUE)

COAST_Study_Data_Clean_Age_Groups_add_dates$covid_average_symptoms_reported <- COAST_Study_Data_Clean_Age_Groups_add_dates$covid_total_symptoms_reported / COAST_Study_Data_Clean_Age_Groups_add_dates$covid_pos_test_num
aggregate(COAST_Study_Data_Clean_Age_Groups_add_dates$covid_average_symptoms_reported,by=list(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group),FUN=mean,na.rm=TRUE)

COAST_Study_Data_Clean_Age_Groups_add_dates$tbi_average_symptoms_reported <- COAST_Study_Data_Clean_Age_Groups_add_dates$tbi_total_symptoms_reported / COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_num
aggregate(COAST_Study_Data_Clean_Age_Groups_add_dates$tbi_average_symptoms_reported,by=list(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group),FUN=mean,na.rm=TRUE)

aggregate(COAST_Study_Data_Clean_Age_Groups_add_dates$covid_pos_test_num,by=list(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group),FUN=mean,na.rm=TRUE)
aggregate(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_tbi_num,by=list(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group),FUN=mean,na.rm=TRUE)
aggregate(COAST_Study_Data_Clean_Age_Groups_add_dates$total_incidences,by=list(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group),FUN=mean,na.rm=TRUE)

#######
# Save outputs (subset of columns needed for downstream: figure_1, efigure_1, etable_*, figure_4)
base_dir   <- if (dir.exists("figure_1")) "figure_1" else "."
output_dir <- file.path(base_dir, "figure_1_output_data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Columns required by downstream scripts (same set as CLEANUP_INPUT_DATA_MINIMAL_COLUMNS.R)
COAST_explicit <- c(
  "participant_id", "qq_group", "recruitment_site", "age_years", "qq_biological_sex",
  "bmi", "qq_education", "race_category", "qq_tbi_num", "covid_pos_test_num", "qq_covid_number",
  "average_years_since_tbi", "average_years_since_covid",
  "recode_sev_score_total_covid", "recode_sev_score_total_tbi", "recode_overall_sev_score",
  "covid_symptom_sev_average", "tbi_symptom_sev_average", "total_symptom_sev_average",
  "recode_freq_score_total_covid", "recode_freq_score_total_tbi", "recode_overall_freq_score",
  "covid_symptom_freq_average", "tbi_symptom_freq_average", "total_symptom_freq_average",
  "qq_phq8_total_score", "qq_phq8_average_score", "qq_eq5d_index_score", "qq_wai_2",
  "qq_fss_average_score", "qq_fss_total_score",
  "qq_gad7_total_score", "chronic_covid", "chronic_tbi", "chronic_acute_overall",
  "recode_tbi_neuro_psych_sev_score", "recode_tbi_neuro_psych_freq_score",
  "recode_covid_neuro_psych_sev_score", "recode_covid_neuro_psych_freq_score",
  "recode_overall_neuro_psych_sev_score", "recode_overall_neuro_psych_freq_score",
  "covid_neuro_psych_symptom_sev_average", "tbi_neuro_psych_symptom_sev_average",
  "total_neuro_psych_symptom_sev_average", "total_neuro_psych_symptom_freq_average",
  "tbi_neuro_psych_symptom_freq_average", "covid_neuro_psych_symptom_freq_average",
  "vib_anti_abeta_1_42_igg_iga", "vib_anti_abeta_1_42_igm",
  "recode_sev_score_total_covid_fatigue", "recode_sev_score_total_covid_headache",
  "recode_sev_score_total_covid_insomnia_sleep_problems", "recode_sev_score_total_covid_drowsiness",
  "recode_sev_score_total_covid_photophobia_phonophobia_se", "recode_sev_score_total_covid_brain_fog",
  "recode_sev_score_total_covid_confusion", "recode_sev_score_total_covid_memory_problems",
  "recode_sev_score_total_covid_difficulty_concentrating", "recode_sev_score_total_covid_difficulty_finding_words",
  "recode_sev_score_total_covid_paresthesia", "recode_sev_score_total_covid_los", "recode_sev_score_total_covid_lot",
  "recode_sev_score_total_covid_dizziness_lightheadedness", "recode_sev_score_total_covid_difficulty_balancing",
  "recode_sev_score_total_tbi_headache", "recode_sev_score_total_tbi_balance_problems",
  "recode_sev_score_total_tbi_dizziness", "recode_sev_score_total_tbi_lightheadedness",
  "recode_sev_score_total_tbi_fatigue", "recode_sev_score_total_tbi_trouble_falling_asleep",
  "recode_sev_score_total_tbi_sleeping_more", "recode_sev_score_total_tbi_drowsiness",
  "recode_sev_score_total_tbi_light_sensitivity", "recode_sev_score_total_tbi_noise_sensitivity",
  "recode_sev_score_total_tbi_irritability", "recode_sev_score_total_tbi_feeling_frustrated_impatient",
  "recode_sev_score_total_tbi_taking_longer_to_think", "recode_sev_score_total_tbi_restlessness",
  "recode_sev_score_total_tbi_sadness", "recode_sev_score_total_tbi_nervousness_anxiousness",
  "recode_sev_score_total_tbi_feeling_more_emotional", "recode_sev_score_total_tbi_feeling_slowed_down",
  "recode_sev_score_total_tbi_in_a_fog", "recode_sev_score_total_tbi_difficulty_concentrating",
  "recode_sev_score_total_tbi_difficulty_remembering", "recode_sev_score_total_tbi_blurred_vision",
  "Age_Group_Tert_Long", "Age_Group_Long", "Age_Group", "dob"
)
d_main <- COAST_Study_Data_Clean_Age_Groups_add_dates
etable2_pattern <- names(d_main)[grepl("^qq_covid_[0-9]+_pos_test_date$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_symptom_onset_date$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_age$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_loc", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_dazed_memory_loss$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_specific_cause$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_penetrating_injury$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_recovery_status$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_med_provider_diagnosis$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_er_doctors_visit$", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_tests___", names(d_main)) |
  grepl("^qq_tbi_[0-9]+_diagnoses___", names(d_main)) |
  grepl("^qq_covid_[0-9]+_pos_test_type$", names(d_main)) |
  grepl("^qq_covid_[0-9]+_med_provider_diagnosis$", names(d_main)) |
  grepl("^qq_covid_[0-9]+_treatments___", names(d_main)) |
  grepl("^qq_covid_[0-9]+_breathing_treatment___", names(d_main)) |
  grepl("^qq_covid_[0-9]+_med_treatments___", names(d_main))]
keep_main <- unique(c(COAST_explicit, etable2_pattern))
keep_main <- keep_main[keep_main %in% names(d_main)]
if (length(keep_main) == 0) keep_main <- COAST_explicit[COAST_explicit %in% names(d_main)]

# Expand to 961 base + 22 symptom columns (983 total) for downstream etables/figures
target_n <- 961L + 22L
if (length(keep_main) < target_n) {
  from_step1 <- setdiff(names(COAST_Study_Data_F100_Clean), keep_main)
  from_step1 <- from_step1[from_step1 %in% names(d_main)]
  add_n <- min(length(from_step1), target_n - length(keep_main))
  if (add_n > 0L) {
    keep_main <- c(keep_main, head(sort(from_step1), add_n))
    keep_main <- keep_main[keep_main %in% names(d_main)]
  }
}

main_csv_path <- normalizePath(file.path(output_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"), mustWork = FALSE)
write.csv(COAST_Study_Data_Clean_Age_Groups_add_dates %>% dplyr::select(dplyr::all_of(keep_main)),
          main_csv_path,
          row.names = FALSE)

# Copy the same CSV to every downstream input folder (absolute paths to avoid working-directory ambiguity)
project_root <- normalizePath(file.path(base_dir, ".."), mustWork = FALSE)
dest_dirs <- c(
  file.path(project_root, "supp_table_1",  "supp_table_1_input_data"),
  file.path(project_root, "supp_table_2",  "supp_table_2_input_data"),
  file.path(project_root, "supp_table_3",  "supp_table_3_input_data"),
  file.path(project_root, "supp_table_4",  "supp_table_4_input_data"),
  file.path(project_root, "supp_table_5",  "supp_table_5_input_data"),
  file.path(project_root, "supp_figure_1", "supp_figure_1_input_data"),
  file.path(project_root, "figure_4",      "figure_4_input_data"),
  file.path(project_root, "supp_figure_5", "supp_figure_5_input_data")
)
for (d in dest_dirs) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ok <- file.copy(main_csv_path, file.path(d, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"), overwrite = TRUE)
  if (!ok) warning("Could not copy CSV to: ", d)
}
message("Copied COAST_Study_Data_Clean_Age_Groups_add_dates.csv to ", length(dest_dirs), " downstream input folders.")

# Only columns needed for COVID positive-test check (no downstream readers)
keep_covid_pos <- c("participant_id", "qq_group",
  paste0("qq_covid_", 1:10, "_test_results"))
keep_covid_pos <- keep_covid_pos[keep_covid_pos %in% names(COAST_Study_Data_Clean_only_covid_pos)]
write.csv(COAST_Study_Data_Clean_only_covid_pos %>% dplyr::select(dplyr::all_of(keep_covid_pos)),
          file.path(output_dir, "COAST_Study_Data_Clean_only_covid_pos.csv"),
          row.names = FALSE)

# Check table: participant_id, qq_group, test results, sum (minimal subset for reproducibility)
check_cols <- c("participant_id", "qq_group", paste0("qq_covid_", 1:10, "_test_results"), "sum_of_pos_tests")
check_cols <- check_cols[check_cols %in% names(COAST_Study_Data_Clean_only_covid_pos_check)]
write.csv(COAST_Study_Data_Clean_only_covid_pos_check %>% dplyr::select(dplyr::all_of(check_cols)),
          file.path(output_dir, "COAST_Study_Data_Clean_only_covid_pos_check.csv"),
          row.names = FALSE)