# =============================================================================
# Supplementary Table 2 – COVID-19 and mTBI incidence-level treatments and tests
# =============================================================================
#
# Description: Reshapes wide COAST data to one row per COVID-19/mTBI incidence.
#   Builds Table 2A (mTBI characteristics) and Table 2B (COVID-19 characteristics).
#   Exported as Word and CSV.
#
# Prerequisites: supp_table_2_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: supp_table_2/supp_table_2A.docx (mTBI), supp_table_2/supp_table_2B.docx (COVID-19);
#   supp_table_2_output_data/supp_table_2A_export.csv, supp_table_2B_export.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) %in% paste0("supp_table_", 1:5)) setwd("..")
input_path <- "supp_table_2/supp_table_2_input_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv"
if (!file.exists(input_path)) stop("Need COAST_Study_Data_Clean_Age_Groups_add_dates.csv in supp_table_2/supp_table_2_input_data/.")
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(input_path, stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(plotrix)
library(gtsummary)
library(flextable)
library(sjstats)
library(stringr)

# -----------------------------------------------------------------------------
# Reshape to long form: one row per COVID-19/mTBI incidence with treatments and tests
# -----------------------------------------------------------------------------

# Convert to long form
all_character_df <- COAST_Study_Data_Clean_Age_Groups_add_dates %>% mutate_all(as.character)
#move columns that we want conserved to the front of the df
all_character_df <- all_character_df %>% relocate(qq_group, Age_Group_Tert_Long, Age_Group_Long, Age_Group, qq_biological_sex, age_years)
# filter down to columns we need
  # matches("^qq_tbi_history___\\d+$"),
  # qq_tbi_num,
  # qq_tbi_multiple,  
wide_to_long_df <- all_character_df %>% dplyr::select(
  participant_id,
  qq_group,
  matches("^qq_covid_\\d+_pos_test_date$"),
  matches("^qq_tbi_\\d+_symptom_onset_date$"),
  matches("^qq_tbi_\\d+_age$"),
  matches("^qq_tbi_\\d+_loc$"),
  matches("^qq_tbi_\\d+_loc_minutes$"),
  matches("^qq_tbi_\\d+_dazed_memory_loss$"),
  matches("^qq_tbi_\\d+_specific_cause$"),
  matches("^qq_tbi_\\d+_penetrating_injury$"),
  matches("^qq_tbi_\\d+_recovery_status$"),
  matches("^qq_tbi_\\d+_med_provider_diagnosis$"),
  matches("^qq_tbi_\\d+_er_doctors_visit$"),
  matches("^qq_tbi_\\d+_tests___\\d+$"),
  matches("^qq_tbi_\\d+_diagnoses___\\d+$"),
  matches("^qq_tbi_\\d+_diagnoses___other$"),
  matches("^qq_covid_\\d+_pos_test_type$"),
  matches("^qq_covid_\\d+_med_provider_diagnosis$"),
  matches("^qq_covid_\\d+_treatments___\\d+$"),
  matches("^qq_covid_\\d+_breathing_treatment___\\d+$"),
  matches("^qq_covid_\\d+_breathing_treatment___other$"),  
  matches("^qq_covid_\\d+_med_treatments___\\d+$"),
  matches("^qq_covid_\\d+_med_treatments___other$"),
)
idx_first <- 3L
idx_last <- ncol(wide_to_long_df)
#long form
long_df_incidences_supp_table_2 <- gather(wide_to_long_df, Label, Value, colnames(wide_to_long_df)[idx_first]:colnames(wide_to_long_df)[idx_last])
# create incidence_num column by grabbing the number from the label
long_df_incidences_supp_table_2$incidence_num <- as.numeric(gsub(".*?([0-9]+).*", "\\1", long_df_incidences_supp_table_2$Label))  
# create incidence_type column by checking what type of incidence a value is from
long_df_incidences_supp_table_2$incidence_type <- ifelse(grepl("tbi", long_df_incidences_supp_table_2$Label), "tbi", "covid")
# remove incidence number from the Labels - make it match
long_df_incidences_supp_table_2 <- long_df_incidences_supp_table_2 %>% mutate(Label = str_replace(Label, paste0(as.character(incidence_num), "_"), ""))
# wide df for incidences
wide_df_incidences_supp_table_2 <- spread(long_df_incidences_supp_table_2, key = Label, value = Value)
# only include rows that actually indicate an incidence
#only include cases where either are not null
wide_df_incidences_supp_table_2_only <- wide_df_incidences_supp_table_2 %>% filter(!is.na(qq_covid_pos_test_date) | !is.na(qq_tbi_symptom_onset_date))
#create incidence date column
wide_df_incidences_supp_table_2_only$incidence_date <- ifelse(!is.na(wide_df_incidences_supp_table_2_only$qq_covid_pos_test_date), wide_df_incidences_supp_table_2_only$qq_covid_pos_test_date, wide_df_incidences_supp_table_2_only$qq_tbi_symptom_onset_date)
wide_df_incidences_supp_table_2_only$incidence_date <- as.Date(wide_df_incidences_supp_table_2_only$incidence_date)
wide_df_incidences_supp_table_2_only <- wide_df_incidences_supp_table_2_only %>% relocate(incidence_date, .after = qq_group)

output_dir <- "supp_table_2/supp_table_2_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Split by incidence type and collapse COVID-19 incidences
# -----------------------------------------------------------------------------
covid_incidences <- wide_df_incidences_supp_table_2_only %>%
  filter(incidence_type == "covid")
# filter out columns with tbi in them
covid_incidences <- covid_incidences %>%
  dplyr::select(-contains("tbi", ignore.case = TRUE))
covid_incidences <- covid_incidences %>%
  dplyr::select(-c("incidence_date", "qq_covid_pos_test_date"))
# create columns for each test type
covid_incidences <- covid_incidences %>%
  mutate(qq_covid_pcr_test_count = ifelse(qq_covid_pos_test_type == "PCR test", 1, 0))
covid_incidences <- covid_incidences %>%
  mutate(qq_covid_antibody_test_count = ifelse(qq_covid_pos_test_type == "Antibody test", 1, 0))
covid_incidences <- covid_incidences %>%
  mutate(qq_covid_rapid_antigen_test_count = ifelse(qq_covid_pos_test_type == "Rapid/antigen test", 1, 0))
covid_incidences <- covid_incidences %>%
  mutate(qq_covid_idk_test_count = ifelse(qq_covid_pos_test_type == "I don't know", 1, 0))
covid_incidences <- covid_incidences %>%
  dplyr::select(-c("incidence_num", "incidence_type", "qq_covid_pos_test_type"))

covid_incidences_numeric <- covid_incidences %>%
  mutate(
    qq_covid_med_provider_diagnosis = case_when(
      is.na(qq_covid_med_provider_diagnosis) ~ 0,  # Convert NA to 0
      qq_covid_med_provider_diagnosis == "No" ~ 0,
      qq_covid_med_provider_diagnosis == "Yes" ~ 1,
      TRUE ~ as.numeric(qq_covid_med_provider_diagnosis)
    )
  )

covid_incidences_numeric <- covid_incidences_numeric %>%
  group_by(participant_id, qq_group) %>%
      mutate_all(~ case_when(
      . == "Unchecked" ~ 0,
      . == "Checked" ~ 1,
      is.na(.) ~ 0,
      TRUE ~ as.numeric(.)
    ))

summarized_covid_incidences_df <- covid_incidences_numeric %>%
  dplyr::group_by(participant_id, qq_group) %>%
  summarise_each(list(sum))

all_cols <- names(summarized_covid_incidences_df)
cols_to_modify <- setdiff(all_cols, c("participant_id", "qq_group"))

transform_checked <- function(x) {
  ifelse(x == 0, "No",
         ifelse(x > 0, "Yes", NA))
}
summarized_covid_incidences_yes_no <- summarized_covid_incidences_df
# Apply the transformation to all selected columns
summarized_covid_incidences_yes_no[cols_to_modify] <- 
  lapply(summarized_covid_incidences_yes_no[cols_to_modify], transform_checked)

# -----------------------------------------------------------------------------
# mTBI incidences
# -----------------------------------------------------------------------------
tbi_incidences <- wide_df_incidences_supp_table_2_only %>%
  filter(incidence_type == "tbi")
# filter out columns with covid in them
tbi_incidences <- tbi_incidences %>%
  dplyr::select(-contains("covid", ignore.case = TRUE))
tbi_incidences <- tbi_incidences %>%
  dplyr::select(-c("incidence_date", "qq_tbi_symptom_onset_date"))
# clean up columns
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_age = as.numeric(str_extract(qq_tbi_age, "\\d+")))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_loc_minutes = as.numeric(qq_tbi_loc_minutes))
# categorical cols
#recovery
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_1_day = ifelse(qq_tbi_recovery_status == "1 day", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_less_than_1_week = ifelse(qq_tbi_recovery_status == "< 1 week", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_1_week_to_1_month = ifelse(qq_tbi_recovery_status == "1 week - 1 month", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_1_month_to_3_months = ifelse(qq_tbi_recovery_status == "1 month - 3 months", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_3_months_to_6_months = ifelse(qq_tbi_recovery_status == "3 months - 6 months", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_6_months_to_1_year = ifelse(qq_tbi_recovery_status == "6 months - 1 year", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_more_than_1_year = ifelse(qq_tbi_recovery_status == ">1 year", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_more_than_5_years = ifelse(qq_tbi_recovery_status == ">5 years", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_still_not_recovered = ifelse(qq_tbi_recovery_status == "I still don't feel recovered", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_recovery_idk = ifelse(qq_tbi_recovery_status == "I don't know", 1, 0))
#cause
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_cause_fall = ifelse(qq_tbi_specific_cause == "Fall (for example, falling from a bike or horse, rollerblading, falling on ice)", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_cause_sports = ifelse(qq_tbi_specific_cause == "Sports", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_cause_fight_trauma = ifelse(qq_tbi_specific_cause == "Fight or trauma (from being hit by someone, or from being shaken violently)", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_cause_car_accident= ifelse(qq_tbi_specific_cause == "Car accident or motor vehicle crash", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_cause_struck_by_object = ifelse(qq_tbi_specific_cause == "Being struck in the head by an object", 1, 0))
tbi_incidences <- tbi_incidences %>%
  mutate(qq_tbi_cause_other= ifelse(qq_tbi_specific_cause == "Other. Please specify:", 1, 0))
# drop not needed cols
tbi_incidences <- tbi_incidences %>%
  dplyr::select(-c("incidence_num", "incidence_type", "qq_tbi_specific_cause", "qq_tbi_recovery_status"))  
# convert things to 1s and 0s
tbi_incidences_numeric <- tbi_incidences %>%
  mutate(
    qq_tbi_dazed_memory_loss = case_when(
      is.na(qq_tbi_dazed_memory_loss) ~ 0,  # Convert NA to 0
      qq_tbi_dazed_memory_loss == "No" ~ 0,
      qq_tbi_dazed_memory_loss == "Yes" ~ 1,
      TRUE ~ as.numeric(qq_tbi_dazed_memory_loss)
    ),
    qq_tbi_er_doctors_visit = case_when(
      is.na(qq_tbi_er_doctors_visit) ~ 0,  # Convert NA to 0
      qq_tbi_er_doctors_visit == "No" ~ 0,
      qq_tbi_er_doctors_visit == "Yes" ~ 1,
      TRUE ~ as.numeric(qq_tbi_er_doctors_visit)
    ),
    qq_tbi_loc = case_when(
      is.na(qq_tbi_loc) ~ 0,  # Convert NA to 0
      qq_tbi_loc == "No" ~ 0,
      qq_tbi_loc == "Yes" ~ 1,
      TRUE ~ as.numeric(qq_tbi_loc)
    ),
    qq_tbi_med_provider_diagnosis = case_when(
      is.na(qq_tbi_med_provider_diagnosis) ~ 0,  # Convert NA to 0
      qq_tbi_med_provider_diagnosis == "No" ~ 0,
      qq_tbi_med_provider_diagnosis == "Yes" ~ 1,
      TRUE ~ as.numeric(qq_tbi_med_provider_diagnosis)
    ),
    qq_tbi_penetrating_injury = case_when(
      is.na(qq_tbi_penetrating_injury) ~ 0,  # Convert NA to 0
      qq_tbi_penetrating_injury == "No" ~ 0,
      qq_tbi_penetrating_injury == "Yes" ~ 1,
      TRUE ~ as.numeric(qq_tbi_penetrating_injury)
    )
  )

tbi_incidences_numeric <- tbi_incidences_numeric %>%
  group_by(participant_id, qq_group) %>%
  mutate_all(
    ~ case_when(
    . == "Unchecked" ~ 0,
    . == "Checked" ~ 1,
    is.na(.) ~ 0,
    TRUE ~ as.numeric(.)
))

summarized_tbi_incidences_df <- tbi_incidences_numeric %>%
  group_by(participant_id, qq_group) %>%
  summarise_each(list(sum))

sum_cols <- tbi_incidences_numeric %>%
  dplyr::select(-qq_tbi_age)
sum_df <- sum_cols %>%
  group_by(participant_id, qq_group) %>%
  summarise_each(list(sum))
mean_cols <- sum_cols <- tbi_incidences_numeric %>%
  dplyr::select(participant_id, qq_group, qq_tbi_age)
mean_df <- mean_cols %>%
  group_by(participant_id, qq_group) %>%
  summarise_each(list(mean))
summarized_tbi_incidences_df <- sum_df %>%
  left_join(mean_df, by = c("participant_id", "qq_group"))

summarized_tbi_incidences_df %>%
  mutate(qq_tbi_dazed_memory_loss = case_when(
    qq_tbi_dazed_memory_loss == 0 ~ "No",
    qq_tbi_dazed_memory_loss > 0 ~ "Yes",
    TRUE ~ as.character(qq_tbi_dazed_memory_loss)
  )) %>% dplyr::select(qq_tbi_dazed_memory_loss)

all_cols <- names(summarized_tbi_incidences_df)
cols_to_modify <- setdiff(all_cols, c("participant_id", "qq_group", "qq_tbi_age", "qq_tbi_loc_minutes"))

transform_checked <- function(x) {
  ifelse(x == 0, "No",
         ifelse(x > 0, "Yes", NA))
}
summarized_tbi_incidences_yes_no <- summarized_tbi_incidences_df
# Apply the transformation to all selected columns
summarized_tbi_incidences_yes_no[cols_to_modify] <- 
  lapply(summarized_tbi_incidences_yes_no[cols_to_modify], transform_checked)


# Table 2A: mTBI characteristics (COVID-19 (-) mTBI (+) vs COVID-19 (+) mTBI (+)); variables and order match paper Table 2
Supplementary_Table_TBI_Data <- summarized_tbi_incidences_yes_no %>%
  ungroup() %>%
  filter(qq_group %in% c("COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)")) %>%
  dplyr::select(qq_group,
                qq_tbi_dazed_memory_loss,
                qq_tbi_age,
                qq_tbi_loc,
                qq_tbi_loc_minutes,
                qq_tbi_er_doctors_visit,
                qq_tbi_med_provider_diagnosis,
                qq_tbi_penetrating_injury,
                qq_tbi_cause_fall,
                qq_tbi_cause_sports,
                qq_tbi_cause_fight_trauma,
                qq_tbi_cause_car_accident,
                qq_tbi_cause_struck_by_object,
                qq_tbi_cause_other,
                qq_tbi_diagnoses___1,
                qq_tbi_diagnoses___2,
                qq_tbi_diagnoses___3,
                qq_tbi_diagnoses___4,
                qq_tbi_diagnoses___5,
                qq_tbi_diagnoses___6,
                qq_tbi_diagnoses___7,
                qq_tbi_diagnoses___8,
                qq_tbi_diagnoses___9,
                qq_tbi_diagnoses___10,
                qq_tbi_diagnoses___11,
                qq_tbi_diagnoses___12,
                qq_tbi_diagnoses___13,
                qq_tbi_diagnoses___14,
                qq_tbi_diagnoses___15,
                qq_tbi_diagnoses___16,
                qq_tbi_tests___1,
                qq_tbi_tests___2,
                qq_tbi_tests___3,
                qq_tbi_tests___4,
                qq_tbi_tests___5,
                qq_tbi_tests___6,
                qq_tbi_recovery_1_day,
                qq_tbi_recovery_less_than_1_week,
                qq_tbi_recovery_1_week_to_1_month,
                qq_tbi_recovery_1_month_to_3_months,
                qq_tbi_recovery_3_months_to_6_months,
                qq_tbi_recovery_6_months_to_1_year,
                qq_tbi_recovery_more_than_1_year,
                qq_tbi_recovery_more_than_5_years,
                qq_tbi_recovery_still_not_recovered,
                qq_tbi_recovery_idk) 


correct_stat_test <- function(data, variable, by, ...) {
  # Check if all values are the same or if there's only one unique value per group
  if (length(unique(data[[variable]])) == 1 || 
      all(tapply(data[[variable]], data[[by]], function(x) length(unique(x))) == 1)) {
    return(tibble(
      test_name = "No test",
      statistic = NA,
      df = NA,
      p.value = NA
    ))
  }
  
  if (is.numeric(data[[variable]])) {
    # For continuous variables, use Kruskal-Wallis or Wilcoxon
    groups <- unique(data[[by]])
    if (length(groups) == 2) {
      # Wilcoxon test for two groups
      test_result <- wilcox.test(data[[variable]] ~ data[[by]])
      test_name <- "Wilcoxon"
    } else {
      # Kruskal-Wallis test for more than two groups
      test_result <- kruskal.test(data[[variable]] ~ data[[by]])
      test_name <- "Kruskal"
    }
    stat <- test_result$statistic
    df <- ifelse(test_name == "Kruskal", test_result$parameter, NA)
    p_value <- test_result$p.value
  } else {
    # For categorical variables, use Chi-squared or Fisher's exact test
    cont_table <- table(data[[by]], data[[variable]])
    expected <- chisq.test(cont_table)$expected
    if (any(expected < 5)) {
      # Fisher's exact test if any expected cell count is less than 5
      test_result <- fisher.test(cont_table, simulate.p.value = TRUE)
      test_name <- "Fisher's"
      stat <- NA
      df <- NA
    } else {
      # Chi-squared test otherwise
      test_result <- chisq.test(cont_table)
      test_name <- "χ²"
      stat <- test_result$statistic
      df <- test_result$parameter
    }
    p_value <- test_result$p.value
  }
  
  # Return a tibble with the results
  tibble(
    test_name = test_name,
    statistic = unname(stat),
    df = unname(df),
    p.value = p_value
  )
}

Supplementary_Table_TBI_Data$qq_group <- factor(Supplementary_Table_TBI_Data$qq_group, levels = c("COVID-19 (-) mTBI (+)", "COVID-19 (+) mTBI (+)"))

Supplementary_Table_TBI <- Supplementary_Table_TBI_Data %>%
  tbl_summary(
    by = qq_group,
    statistic = list(all_continuous() ~ "{mean} ({std.error})", all_categorical() ~ "{n} ({p}%)"),
    missing = "no",
    digits = list(all_categorical() ~ c(0, 1), all_continuous() ~ c(2, 2)),
    label = list(qq_tbi_dazed_memory_loss ~ "Dazed/Gap in memory",
                 qq_tbi_age ~ "Average age of injury",
                 qq_tbi_loc ~ "Loss of consciousness",
                 qq_tbi_loc_minutes ~ "Loss of consciousness (minutes)",
                 qq_tbi_er_doctors_visit ~ "ER visit",
                 qq_tbi_med_provider_diagnosis ~ "Injury diagnosed by medical provider",
                 qq_tbi_penetrating_injury ~ "Penetrating injury",
                 qq_tbi_cause_fall ~ "Fall",
                 qq_tbi_cause_sports ~ "Sports",
                 qq_tbi_cause_fight_trauma ~ "Fight/Trauma",
                 qq_tbi_cause_car_accident ~ "Motor vehicle accident",
                 qq_tbi_cause_struck_by_object ~ "Being struck in the head by an object",
                 qq_tbi_cause_other ~ "Other cause",
                 qq_tbi_diagnoses___1 ~ "Subdural hematoma",
                 qq_tbi_diagnoses___2 ~ "Epidural hematoma",
                 qq_tbi_diagnoses___3 ~ "Skull fracture",
                 qq_tbi_diagnoses___4 ~ "Hygroma",
                 qq_tbi_diagnoses___5 ~ "Intraventricular hemorrhage",
                 qq_tbi_diagnoses___6 ~ "Subarachnoid hemorrhage",
                 qq_tbi_diagnoses___7 ~ "Brain Contusion",
                 qq_tbi_diagnoses___8 ~ "Axonal Injury",
                 qq_tbi_diagnoses___9 ~ "Increased intracranial pressure (ICP)",
                 qq_tbi_diagnoses___10 ~ "Cerebral Edema",
                 qq_tbi_diagnoses___11 ~ "Cerebral herniation",
                 qq_tbi_diagnoses___12 ~ "Cerebral ischemia",
                 qq_tbi_diagnoses___13 ~ "Vertebral artery occlusion",
                 qq_tbi_diagnoses___14 ~ "Other vascular injury",
                 qq_tbi_diagnoses___15 ~ "Other",
                 qq_tbi_diagnoses___16 ~ "None of the above",
                 qq_tbi_tests___1 ~ "Computerized tomography (CT) scan",
                 qq_tbi_tests___2 ~ "Magnetic resonance imaging (MRI) scan",
                 qq_tbi_tests___3 ~ "Positron emission tomography (PET) scan",
                 qq_tbi_tests___4 ~ "Single-photon emission computed tomography (SPECT) scan",
                 qq_tbi_tests___5 ~ "Angiography",
                 qq_tbi_tests___6 ~ "None of the above",
                 qq_tbi_recovery_1_day ~ "1 day",
                 qq_tbi_recovery_less_than_1_week ~ "<1 week",
                 qq_tbi_recovery_1_week_to_1_month ~ "1 week - 1 month",
                 qq_tbi_recovery_1_month_to_3_months ~ "1 month - 3 months",
                 qq_tbi_recovery_3_months_to_6_months ~ "3 months - 6 months",
                 qq_tbi_recovery_6_months_to_1_year ~ "6 months - 1 year",
                 qq_tbi_recovery_more_than_1_year ~ ">1 year",
                 qq_tbi_recovery_more_than_5_years ~ ">5 years",
                 qq_tbi_recovery_still_not_recovered ~ "Still don't feel recovered at time of study visit",
                 qq_tbi_recovery_idk ~ "Don't know")) %>%
  add_stat(fns = everything() ~ correct_stat_test) %>%
  modify_header(
    test_name = "**Test**",
    statistic = "**Test Statistic**",
    df = "**Degrees of Freedom**",
    p.value = "**p-value**"
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05) %>%
  modify_footnote(
    update = starts_with("stat_") ~
      "Mean (SEM) for continuous variables; n (%) for categorical variables"
  )

Supplementary_Table_TBI_Final <- as_flex_table(Supplementary_Table_TBI)
Supplementary_Table_TBI_Final

save_as_docx(Supplementary_Table_TBI_Final, path = "supp_table_2/supp_table_2A.docx")
write.csv(Supplementary_Table_TBI$table_body, file.path(output_dir, "supp_table_2A_export.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# Table 2B: COVID-19 incidence characteristics (COVID-19 (+) mTBI (-) vs COVID-19 (+) mTBI (+))
# -----------------------------------------------------------------------------
Supplementary_Table_COVID <- summarized_covid_incidences_yes_no %>%
  ungroup() %>%
  filter(qq_group %in% c("COVID-19 (+) mTBI (-)", "COVID-19 (+) mTBI (+)")) %>%
  dplyr::select(qq_group,          
                                                                      qq_covid_pcr_test_count,           
                                                                      qq_covid_antibody_test_count,  
                                                                      qq_covid_rapid_antigen_test_count, 
                                                                      qq_covid_idk_test_count, 
                                                                      qq_covid_med_provider_diagnosis,
                                                                      qq_covid_treatments___2,               
                                                                      qq_covid_treatments___3,              
                                                                      # qq_covid_breathing_treatment___1,    
                                                                      qq_covid_breathing_treatment___2,     
                                                                      qq_covid_breathing_treatment___3,     
                                                                      qq_covid_breathing_treatment___4,    
                                                                      qq_covid_breathing_treatment___5,     
                                                                      qq_covid_med_treatments___1, 
                                                                      qq_covid_med_treatments___2,           
                                                                      qq_covid_med_treatments___3,         
                                                                      qq_covid_med_treatments___4,          
                                                                      qq_covid_med_treatments___5,          
                                                                      qq_covid_med_treatments___6,         
                                                                      qq_covid_med_treatments___7,          
                                                                      qq_covid_med_treatments___8,          
                                                                      qq_covid_med_treatments___9,   
                                                                      qq_covid_med_treatments___10,         
                                                                      qq_covid_med_treatments___11,         
                                                                      qq_covid_med_treatments___12,        
                                                                      qq_covid_med_treatments___13) 


Supplementary_Table_COVID_clean <- Supplementary_Table_COVID

Supplementary_Table_COVID_clean$hospitalized <- ifelse(Supplementary_Table_COVID$qq_covid_breathing_treatment___3 == "Yes" | 
                                                       Supplementary_Table_COVID$qq_covid_breathing_treatment___4 == "Yes", 
                        "Yes", "No")

Supplementary_Table_COVID_clean$qq_group <- factor(Supplementary_Table_COVID_clean$qq_group, levels = c("COVID-19 (+) mTBI (-)", "COVID-19 (+) mTBI (+)"))
Supplementary_Table_COVID <- Supplementary_Table_COVID_clean %>%
  tbl_summary(
    by = qq_group,
    statistic = list(all_continuous() ~ "{mean} ({std.error})", all_categorical() ~ "{n} ({p}%)"),
    missing = "no",
    digits = list(all_categorical() ~ c(0, 2), all_continuous() ~ c(2, 2)),
    label = list(qq_covid_pcr_test_count ~ "PCR test",           
                 qq_covid_antibody_test_count ~ "Antibody test",  
                 qq_covid_rapid_antigen_test_count ~ "Rapid/Antigen test", 
                 qq_covid_idk_test_count ~ "Test type unknown", 
                 qq_covid_med_provider_diagnosis ~ "Infection diagnosed by medical provider", 
                 # qq_covid_breathing_treatment___1 ~ "I did not receive breathing treatment",    
                 qq_covid_breathing_treatment___2 ~ "Oxygen (through an oxygen mask or tube under my nose, no pressure applied)",     
                 qq_covid_breathing_treatment___3 ~ "Oxygen (through an oxygen mask, which pushes oxygen into your lungs)",     
                 qq_covid_breathing_treatment___4 ~ "A breathing machine (ventilator) with a tube down my throat",    
                 qq_covid_breathing_treatment___5 ~ "Other breathing treatment",     
                 qq_covid_med_treatments___1 ~ "Oral corticosteroids (eg. Prednisone, Dexamethasone, Methylprednisolone, Hydrocortisone)", 
                 qq_covid_med_treatments___2 ~ "Inhaled corticosteroids (eg. fluticasone (Flovent), beclomethasone (QVar), etc.)",           
                 qq_covid_med_treatments___3 ~ "Anti-TNF medications (infliximab, adalimumab, certolizumab, golimumab, etanercept, others)",         
                 qq_covid_med_treatments___4 ~ "IL-6 pathway inhibitors (sarilumab,tocilizumab, siltuximab, others)",          
                 qq_covid_med_treatments___5 ~ "Non-steroidal anti-inflammatory agents (NSAIDS) with or without a prescription: (eg. ibuprofen (Motrin, Advil), naproxen (Naprosyn, Aleve, Anaprox, Naprelan), diclofenac (Cambia, Cataflam, Voltaren, Zipsor), indomethacin (Indocin), diflunisal, etodolac, ketoprofen, ketorolac, nambumetone, oxaprozin (Daypro), piroxicam (Feldene), salsalate (Disalate), sulidnac, tolmetin, celecoxib (Celebrex)",          
                 qq_covid_med_treatments___6 ~ "Antiviral drugs (Remdesivir, Paxlovid)",         
                 qq_covid_med_treatments___7 ~ "Antiparasitic drugs (Ivermectin)",          
                 qq_covid_med_treatments___8 ~ "Antibiotic drugs (azithromycin, doxycycline, clarithromycin, ceftriaxone, erythromycin, amoxicillin, ampicillin, gentamicin, benzylpenicillin, piperacillin/tazobactam, ciprofloxacin, ceftazidime, cefepime, vancomycin, meropenem, cefuroxime, etc.)",          
                 qq_covid_med_treatments___9 ~ "Hydroxychloroquine or chloroquine",   
                 qq_covid_med_treatments___10 ~ "Kinase inhibitors (Acalabrutinib (Calquence); Baricitinib (Olumiant); Ruxolitinib (Jakafi); Tofacitinib (Xeljanz)",         
                 qq_covid_med_treatments___11 ~ "Monoclonal antibody (MABs) treatment",         
                 qq_covid_med_treatments___12 ~ "Other",        
                 qq_covid_med_treatments___13 ~ "None of the above", 
                 qq_covid_treatments___2 ~ "I recovered at home",               
                 qq_covid_treatments___3 ~ "I spoke with a healthcare professional and wasn’t admitted to the hospital",              
                 hospitalized ~ "I was admitted to the hospital for at least one night")) %>%
  add_stat(fns = everything() ~ correct_stat_test) %>%
  modify_header(
    test_name = "**Test**",
    statistic = "**Test Statistic**",
    df = "**Degrees of Freedom**",
    p.value = "**p-value**"
  ) %>%
  bold_p(t = 0.05) %>%
  modify_header(label = "** **") %>%
  modify_footnote(
    update = starts_with("stat_") ~
      "Mean (SEM) for continuous variables; n (%) for categorical variables") %>%
  bold_labels()

Supplementary_Table_COVID_Final <- as_flex_table(Supplementary_Table_COVID)
Supplementary_Table_COVID_Final

save_as_docx(Supplementary_Table_COVID_Final, path = "supp_table_2/supp_table_2B.docx")
write.csv(Supplementary_Table_COVID$table_body, file.path(output_dir, "supp_table_2B_export.csv"), row.names = FALSE)

