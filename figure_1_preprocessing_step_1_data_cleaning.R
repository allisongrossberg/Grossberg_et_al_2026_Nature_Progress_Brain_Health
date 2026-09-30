# =============================================================================
# COAST Study Year 1 – Qualtrics data cleaning (first 100 participants)
# =============================================================================
#
# Description: Cleans and recodes COAST Year 1 Qualtrics export for analysis:
#   demographics (height, age, race), BMI, PHQ-8, GAD-7, FSS, PSS, Neuro-QoL,
#   and APOE result. Excludes test records (participant_id containing "DBSV1").
#
# Input:  figure_1_input_data/COAST_Data_Y1_8_1_24.csv
#         figure_1_input_data/COAST_Headers_Y1_8_1_24.csv
# Output: figure_1_output_data/COAST_Study_Data_F100_Clean.csv
#         (Full cleaned dataset required by step 2; step 2 then writes subset CSVs for downstream.)
#
# Usage:  Set working directory to figure_1, then run this script.
#
# Packages: tidyr, dplyr, lubridate, readxl, readr, stringr (tidyverse)
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(lubridate)
library(readxl)
library(readr)
library(stringr)

# -----------------------------------------------------------------------------
# Paths (run from figure_1 or project root)
# -----------------------------------------------------------------------------
base_dir   <- if (dir.exists("figure_1")) "figure_1" else "."
input_dir  <- file.path(base_dir, "figure_1_input_data")
output_dir <- file.path(base_dir, "figure_1_output_data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

path_data    <- file.path(input_dir, "COAST_Data_Y1_8_1_24.csv")
path_headers <- file.path(input_dir, "COAST_Headers_Y1_8_1_24.csv")
path_output  <- file.path(output_dir, "COAST_Study_Data_F100_Clean.csv")

# -----------------------------------------------------------------------------
# Load data and apply clean column names
# -----------------------------------------------------------------------------
COAST_Data_Y1_8_1_24 <- read_csv(path_data, show_col_types = FALSE)
headers_df           <- read_csv(path_headers, show_col_types = FALSE)
colnames(COAST_Data_Y1_8_1_24) <- colnames(headers_df)

COAST_Study_Data_F100 <- COAST_Data_Y1_8_1_24

# -----------------------------------------------------------------------------
# Filter to analysis sample (exclude test records)
# -----------------------------------------------------------------------------
nrow(COAST_Study_Data_F100)
COAST_Study_Data_F100_Clean <- COAST_Study_Data_F100 %>% filter(!grepl('DBSV1', participant_id))
nrow(COAST_Study_Data_F100_Clean)

# -----------------------------------------------------------------------------
# Clean height column (original format: feet'inches)
# ----------------------------------------------------------------------------- 

COAST_Study_Data_F100_Clean$height_feet <- as.numeric(str_split_fixed(COAST_Study_Data_F100_Clean$qq_height, "'", 3)[,1])
COAST_Study_Data_F100_Clean$qq_height <- as.numeric(str_split_fixed(COAST_Study_Data_F100_Clean$qq_height, "'", 3)[,2])
COAST_Study_Data_F100_Clean$height_inches_final <- COAST_Study_Data_F100_Clean$qq_height + (COAST_Study_Data_F100_Clean$height_feet * 12) 

# -----------------------------------------------------------------------------
# Age at survey start
# -----------------------------------------------------------------------------
# age_years arrives pre-computed in the raw export: exact date of birth is a
# direct identifier and is not included in this deposit (age_years was
# computed once from DOB + survey start date, then DOB was removed).

COAST_Study_Data_F100_Clean$age_years <- as.numeric(COAST_Study_Data_F100_Clean$age_years)

# qq_start_date contains a date and time; split into separate columns.
COAST_Study_Data_F100_Clean <- COAST_Study_Data_F100_Clean %>% separate(qq_start_date, c('qq_start_date', 'qq_start_time'), sep = " ")
COAST_Study_Data_F100_Clean$qq_start_date <- as.Date(COAST_Study_Data_F100_Clean$qq_start_date, format = "%Y-%m-%d")

# -----------------------------------------------------------------------------
# Race category (from checkbox fields)
# -----------------------------------------------------------------------------

COAST_Study_Data_F100_Clean$race_category <- NA
COAST_Study_Data_F100_Clean$race_category <- ifelse(COAST_Study_Data_F100_Clean$qq_race___1 == "Checked", "White", NA)
COAST_Study_Data_F100_Clean$race_category <- ifelse(COAST_Study_Data_F100_Clean$qq_race___2 == "Checked" & is.na(COAST_Study_Data_F100_Clean$race_category), "Black, African American", ifelse(COAST_Study_Data_F100_Clean$qq_race___2 == "Checked" & is.character(COAST_Study_Data_F100_Clean$race_category), "Two or more races", COAST_Study_Data_F100_Clean$race_category))
COAST_Study_Data_F100_Clean$race_category <- ifelse(COAST_Study_Data_F100_Clean$qq_race___3 == "Checked" & is.na(COAST_Study_Data_F100_Clean$race_category), "Native or indigenous to U.S. lands", ifelse(COAST_Study_Data_F100_Clean$qq_race___2 == "Checked" & is.character(COAST_Study_Data_F100_Clean$race_category), "Two or more races", COAST_Study_Data_F100_Clean$race_category))
COAST_Study_Data_F100_Clean$race_category <- ifelse(COAST_Study_Data_F100_Clean$qq_race___4 == "Checked" & is.na(COAST_Study_Data_F100_Clean$race_category), "Asian", ifelse(COAST_Study_Data_F100_Clean$qq_race___2 == "Checked" & is.character(COAST_Study_Data_F100_Clean$race_category), "Two or more races", COAST_Study_Data_F100_Clean$race_category))
COAST_Study_Data_F100_Clean$race_category <- ifelse(grepl("Hispanic", COAST_Study_Data_F100_Clean$qq_race_other, fixed = TRUE), "Hispanic or Latino", COAST_Study_Data_F100_Clean$race_category)
COAST_Study_Data_F100_Clean$race_category <- ifelse(grepl("Mixed", COAST_Study_Data_F100_Clean$qq_race_other, fixed = TRUE), "Two or more races", COAST_Study_Data_F100_Clean$race_category)
COAST_Study_Data_F100_Clean$race_category <- ifelse(is.na(COAST_Study_Data_F100_Clean$race_category), "Declined to respond", COAST_Study_Data_F100_Clean$race_category)

# -----------------------------------------------------------------------------
# Calculate BMI
# ----------------------------------------------------------------------------- 

COAST_Study_Data_F100_Clean <- COAST_Study_Data_F100_Clean %>% mutate(bmi = (703 * qq_weight_lbs) / (height_inches_final^2))

# -----------------------------------------------------------------------------
# Clean PHQ-8 (recode to numeric)
# ----------------------------------------------------------------------------- 

phq8_recode <- grep("_phq8_[0-9]+", names(COAST_Study_Data_F100_Clean), value = TRUE)
phq8_recode

for (col_name in phq8_recode) {
  COAST_Study_Data_F100_Clean[paste0(col_name, "_recode")] <- ifelse(COAST_Study_Data_F100_Clean[col_name] == "Not at all", 0, 
                                                                    ifelse(COAST_Study_Data_F100_Clean[col_name] == "Several days", 1, 
                                                                           ifelse(COAST_Study_Data_F100_Clean[col_name] == "More than half the days", 2, 
                                                                                  ifelse(COAST_Study_Data_F100_Clean[col_name] == "Nearly every day", 3, 
                                                                                         ifelse(COAST_Study_Data_F100_Clean[col_name] %in% c("Not at all","Several days", "More than half the days", "Nearly every day"), COAST_Study_Data_F100_Clean[col_name],NA)))))
}

# -----------------------------------------------------------------------------
# Clean GAD-7 (recode to numeric)
# -----------------------------------------------------------------------------

gad7recode <- grep("_gad7_[0-9]+", names(COAST_Study_Data_F100_Clean), value = TRUE)
gad7recode

for (col_name in gad7recode) {
  COAST_Study_Data_F100_Clean[paste0(col_name, "_recode")] <- ifelse(COAST_Study_Data_F100_Clean[col_name] == "Not at all", 0, 
                                                                    ifelse(COAST_Study_Data_F100_Clean[col_name] == "Several days", 1, 
                                                                           ifelse(COAST_Study_Data_F100_Clean[col_name] == "More than half the days", 2, 
                                                                                  ifelse(COAST_Study_Data_F100_Clean[col_name] == "Nearly every day", 3, 
                                                                                         ifelse(COAST_Study_Data_F100_Clean[col_name] %in% c("Not at all","Several days", "More than half the days", "Nearly every day"), COAST_Study_Data_F100_Clean[col_name],NA)))))
}

# -----------------------------------------------------------------------------
# Clean FSS (recode to numeric)
# -----------------------------------------------------------------------------

fssrecode <- grep("_fss_[0-9]+", names(COAST_Study_Data_F100_Clean), value = TRUE)
fssrecode

for (col_name in fssrecode) {
  COAST_Study_Data_F100_Clean[paste0(col_name, "_recode")] <- ifelse(COAST_Study_Data_F100_Clean[col_name] == "Strongly Disagree", 1, 
                                                                    ifelse(COAST_Study_Data_F100_Clean[col_name] == "Disagree", 2, 
                                                                           ifelse(COAST_Study_Data_F100_Clean[col_name] == "Somewhat Disagree", 3, 
                                                                                  ifelse(COAST_Study_Data_F100_Clean[col_name] == "Neutral", 4, 
                                                                                         ifelse(COAST_Study_Data_F100_Clean[col_name] == "Somewhat Agree", 5, 
                                                                                                ifelse(COAST_Study_Data_F100_Clean[col_name] == "Agree", 6, 
                                                                                                       ifelse(COAST_Study_Data_F100_Clean[col_name] == "Strongly Agree", 7, 
                                                                                                              ifelse(COAST_Study_Data_F100_Clean[col_name] %in% c("Strongly Disagree", "Disagree", "Somewhat Disagree", "Neutral","Somewhat Agree", "Agree", "Strongly Agree"), 
                                                                                                                     COAST_Study_Data_F100_Clean[col_name],NA))))))))
}

# -----------------------------------------------------------------------------
# Clean PSS (recode to numeric)
# -----------------------------------------------------------------------------

COAST_Study_Data_F100_Clean$qq_pss_1_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_1 == "Never", 0, ifelse(COAST_Study_Data_F100_Clean$qq_pss_1== "Almost Never", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_1== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_1== "Fairly Often", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_1== "Very Often", 4, NA))))) 
COAST_Study_Data_F100_Clean$qq_pss_2_recode<- ifelse(COAST_Study_Data_F100_Clean$qq_pss_2 == "Never", 0, ifelse(COAST_Study_Data_F100_Clean$qq_pss_2== "Almost Never", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_2== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_2== "Fairly Often", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_2== "Very Often", 4, NA)))))
COAST_Study_Data_F100_Clean$qq_pss_3_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_3 == "Never", 0, ifelse(COAST_Study_Data_F100_Clean$qq_pss_3== "Almost Never", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_3== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_3== "Fairly Often", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_3== "Very Often", 4, NA)))))
COAST_Study_Data_F100_Clean$qq_pss_4_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_4 == "Never", 4, ifelse(COAST_Study_Data_F100_Clean$qq_pss_4== "Almost Never", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_4== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_4== "Fairly Often", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_4== "Very Often", 0, NA)))))
COAST_Study_Data_F100_Clean$qq_pss_5_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_5 == "Never", 4, ifelse(COAST_Study_Data_F100_Clean$qq_pss_5== "Almost Never", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_5== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_5== "Fairly Often", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_5== "Very Often", 0, NA)))))
COAST_Study_Data_F100_Clean$qq_pss_6_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_6 == "Never", 0, ifelse(COAST_Study_Data_F100_Clean$qq_pss_6== "Almost Never", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_6== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_6== "Fairly Often", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_6== "Very Often", 4, NA))))) 
COAST_Study_Data_F100_Clean$qq_pss_7_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_7 == "Never", 4, ifelse(COAST_Study_Data_F100_Clean$qq_pss_7== "Almost Never", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_7== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_7== "Fairly Often", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_7== "Very Often", 0, NA)))))
COAST_Study_Data_F100_Clean$qq_pss_8_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_8 == "Never", 4, ifelse(COAST_Study_Data_F100_Clean$qq_pss_8== "Almost Never", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_8== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_8== "Fairly Often", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_8== "Very Often", 0, NA))))) 
COAST_Study_Data_F100_Clean$qq_pss_9_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_9 == "Never", 0, ifelse(COAST_Study_Data_F100_Clean$qq_pss_9== "Almost Never", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_9== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_9== "Fairly Often", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_9== "Very Often", 4, NA))))) 
COAST_Study_Data_F100_Clean$qq_pss_10_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_pss_10 == "Never", 0, ifelse(COAST_Study_Data_F100_Clean$qq_pss_10== "Almost Never", 1, ifelse(COAST_Study_Data_F100_Clean$qq_pss_10== "Sometimes", 2, ifelse(COAST_Study_Data_F100_Clean$qq_pss_10== "Fairly Often", 3, ifelse(COAST_Study_Data_F100_Clean$qq_pss_10== "Very Often", 4, NA))))) 

# -----------------------------------------------------------------------------
# Clean Neuro-QoL (recode to numeric)
# ----------------------------------------------------------------------------- 

COAST_Study_Data_F100_Clean$qq_neuroqol_1_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_1 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_1== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_1== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_1== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_1== "Cannot do", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_2_recode<- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_2 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_2== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_2== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_2== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_2== "Cannot do", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_3_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_3 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_3== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_3== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_3== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_3== "Cannot do", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_4_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_4 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_4== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_4== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_4== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_4== "Cannot do", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_5_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_5 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_5== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_5== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_5== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_5== "Cannot do", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_6_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_6 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_6== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_6== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_6== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_6== "Cannot do", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_7_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_7 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_7== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_7== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_7== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_7== "Cannot do", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_8_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_8 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_8== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_8== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_8== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_8== "Cannot do", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_9_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_9 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_9== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_9== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_9== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_9== "Cannot do", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_10_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_10 == "None", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_10== "A little", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_10== "Somewhat", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_10== "A lot", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_10== "Cannot do", 1, NA))))) 

COAST_Study_Data_F100_Clean$qq_neuroqol_11_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_11 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_11== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_11== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_11== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_11== "Very often (several times a day)", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_12_recode<- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_12 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_12== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_12== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_12== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_12== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_13_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_13 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_13== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_13== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_13== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_13== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_14_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_14 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_14== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_14== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_14== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_14== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_15_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_15 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_15== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_15== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_15== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_15== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_16_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_16 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_16== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_16== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_16== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_16== "Very often (several times a day)", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_17_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_17 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_17== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_17== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_17== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_17== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_18_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_18 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_18== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_18== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_18== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_18== "Very often (several times a day)", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_19_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_19 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_19== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_19== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_19== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_19== "Very often (several times a day)", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_20_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_20 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_20== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_20== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_20== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_20== "Very often (several times a day)", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_21_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_21 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_21== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_21== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_21== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_21== "Very often (several times a day)", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_22_recode<- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_22 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_22== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_22== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_22== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_22== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_23_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_23 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_23== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_23== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_23== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_23== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_24_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_24 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_24== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_24== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_24== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_24== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_25_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_25 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_25== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_25== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_25== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_25== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_26_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_26 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_26== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_26== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_26== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_26== "Very often (several times a day)", 1, NA))))) 
COAST_Study_Data_F100_Clean$qq_neuroqol_27_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_27 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_27== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_27== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_27== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_27== "Very often (several times a day)", 1, NA)))))
COAST_Study_Data_F100_Clean$qq_neuroqol_28_recode <- ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_28 == "Never", 5, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_28== "Rarely (once)", 4, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_28== "Sometimes (2-3 times)", 3, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_28== "Often (once a day)", 2, ifelse(COAST_Study_Data_F100_Clean$qq_neuroqol_28== "Very often (several times a day)", 1, NA))))) 

# -----------------------------------------------------------------------------
# APOE result string normalization
# -----------------------------------------------------------------------------
COAST_Study_Data_F100_Clean$vib_apoe_result <- str_replace(COAST_Study_Data_F100_Clean$vib_apoe_result, "E2/3", "E2/E3")

write.csv(COAST_Study_Data_F100_Clean, path_output, row.names = FALSE)

# -----------------------------------------------------------------------------
# Data quality note (symptom dates)
# -----------------------------------------------------------------------------
# For some participants who took the Year 1 survey in 2024, COVID-19 symptom
# onset and resolution dates could not be reported as 2024 due to a Qualtrics
# constraint; many reported 2023 instead. No correction is applied, as true
# intent (2023 vs 2024) cannot be determined.
