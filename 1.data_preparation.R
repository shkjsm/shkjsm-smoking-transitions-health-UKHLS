################################################################################
# Dissertation Project: The Effect of Smokers transitioning to E-cigarette in
#                       Physical and Mental health : An Emulated Trial using Longitudinal Data
#
# SCRIPT 1: DATA PREPARATION
#
# What this script does:
# 1. Loads the raw Stata dataset.
# 2. Recodes all necessary variables into a consistent format.
# 3. Calculates real household income adjusted for inflation.
# 4. Creates and saves a master "complete cases" dataset for analysis.
#
################################################################################
#===============================================================================
# 1. SETUP
#===============================================================================
library(tidyverse) # For data manipulation and plotting
library(readstata13) # To read Stata v13 files

# --- Define paths ---
DATA_PATH <- "C:/Users/shaik/OneDrive - University of Glasgow/MPH/Research Project/Project/Source data"
OUTPUT_DIR <- "./output_data/"

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE)

#===============================================================================
# 2. LOAD AND CLEAN DATA
#===============================================================================

# --- Load the raw Stata dataset ---
ukhls_raw <- read.dta13(paste0(DATA_PATH, "/jasim_linked.dta"))

# --- Master cleaning and recoding pipe ---
ukhls_clean <- ukhls_raw %>%
  mutate(
    # Recode job status
    jbstat = case_when(
      jbstat %in% c("1. Self employed", "2. Paid employment(ft/pt)",  "9. Govt training scheme", 
                    "10. Unpaid, family business", "11. On apprenticeship", "12. On furlough", 
                    "13. Temporarily laid off/short term working",
                    "14. On shared parental leave", "15. On adoption leave") ~ "Employed",
      
      
      jbstat %in% c("3. Unemployed") ~ "Unemployed",
      
      jbstat %in% c("4. Retired", "5. On maternity leave", "6. Family care or home", "7. Full-time student",
                    "8. LT sick or disabled", "97. Doing something else") ~ "Economically inactive",
      TRUE ~ NA_character_
    ),
    # Recode long-term health condition
    health = case_when(health == "1. Yes" ~ "Yes", health == "2. No" ~ "No", TRUE ~ NA_character_),
    
    # Recode GP visits
    hl2gp = case_when(
      hl2gp == "0. None" ~ "Zero",
      hl2gp %in% c("1. One or two", "2. Three to five") ~ "One to five",
      hl2gp == "3. Six to ten" ~ "Six to ten",
      hl2gp == "4. More than ten" ~ "More than ten",
      TRUE ~ NA_character_
    ),
    
    # Recode smoking status
    smoker = case_when(smoker == "1. Yes" ~ "Smoker", smoker == "2. No" ~ "Non-smoker", TRUE ~ NA_character_),
    
    # Recode e-cigarette use
    ecigs_use = coalesce(ecigs, ecigs1),
    ecigs_use = case_when(
      ecigs_use %in% c("1. I have never used e-cigarettes", "2. I have only tried using e-cigarettes once or twice",
                       "3. I used e-cigarettes regularly in the past, but I never use them now") ~ "Non-user",
      
      ecigs_use %in% c("5. I use e-cigarettes at least once a month, but less than once a month",
                       "4. I sometimes use e-cigarettes but less than once a month",
                       "6. I use e-cigarettes at least once a week") ~ "Current user",
      TRUE ~ NA_character_
    ),
    
    # Recode sex
    sex_dv = case_when(sex_dv == "1. Male" ~ "Male", sex_dv == "2. Female" ~ "Female", TRUE ~ NA_character_),
    
    # Clean age and create squared term
    age_dv = na_if(age_dv, -9),
    age_squared = age_dv^2,
    
    # Recode ethnicity
    ethn_dv = case_when(
      ethn_dv %in% c("1. british/english/scottish/welsh/northern irish", "2. irish", "3. gypsy or irish traveller", "4. any other white background") ~ "White",
      ethn_dv %in% c("14. caribbean", "15. african", "16. any other black background") ~ "Black",
      ethn_dv %in% c("5. white and black caribbean", "6. white and black african", "7. white and asian", "8. any other mixed background") ~ "Mixed",
      ethn_dv %in% c("9. indian", "10. pakistani", "11. bangladeshi", "12. chinese", "13. any other asian background") ~ "Asian",
      ethn_dv %in% c("17. arab", "97. any other ethnic group") ~ "Other",
      TRUE ~ NA_character_
    ),
    
    # Recode highest qualification
    hiqual_dv = case_when(
      hiqual_dv %in% c("1. Degree", "2. Other higher degree") ~ "Degree or higher",
      hiqual_dv == "3. A-level etc" ~ "A-level or equivalent",
      hiqual_dv == "4. GCSE etc" ~ "GCSE or equivalent",
      hiqual_dv == "5. Other qualification" ~ "Other qualification",
      hiqual_dv == "9. No qualification" ~ "No qualification",
      TRUE ~ NA_character_
    ),
    
    # Recode region
    gor_dv = case_when(
      gor_dv == "1. North East" ~ "North East",
      gor_dv == "2. North West" ~ "North West",
      gor_dv == "3. Yorkshire and the Humber" ~ "Yorkshire and the Humber",
      gor_dv == "4. East Midlands" ~ "East Midlands",
      gor_dv == "5. West Midlands" ~ "West Midlands",
      gor_dv == "6. East of England" ~ "East of England",
      gor_dv == "7. London" ~ "London",
      gor_dv == "8. South East" ~ "South East",
      gor_dv == "9. South West" ~ "South West",
      gor_dv == "10. Wales" ~ "Wales",
      gor_dv == "11. Scotland" ~ "Scotland",
      gor_dv == "12. Northern Ireland" ~ "Northern Ireland",
      TRUE ~ NA_character_
    ),
    
    # Clean health scores and income (replace negative codes with NA)
    sf12pcs_dv = ifelse(sf12pcs_dv < 0, NA, sf12pcs_dv),
    sf12mcs_dv = ifelse(sf12mcs_dv < 0, NA, sf12mcs_dv),
    fihhmnnet1_dv = ifelse(fihhmnnet1_dv < 0, NA, fihhmnnet1_dv),
    ieqmoecd_dv = ifelse(ieqmoecd_dv < 0, NA, ieqmoecd_dv)
  )
# --- Calculate inflation-adjusted real income ---
cpi_data <- tibble(wave = 7:14, cpi = c(103.6, 106.0, 107.8, 108.9, 111.6, 120.5, 128.6, 132.9))

ukhls_final <- ukhls_clean %>%
  mutate(equivalised_income = fihhmnnet1_dv / ieqmoecd_dv) %>%
  left_join(cpi_data, by = "wave") %>%
  mutate(
    real_income = equivalised_income * (100 / cpi),
    log_real_income = ifelse(real_income > 0, log(real_income), NA)
  )

# --- Create final complete case dataset ---
ukhls_complete_cases <- ukhls_final %>%
  select(
    pidp, wave, jbstat, health, hl2gp, smoker, ecigs_use, sex_dv, age_dv,
    age_squared, ethn_dv, gor_dv, hiqual_dv, sf12pcs_dv, sf12mcs_dv,
    nkids015, log_real_income
  ) %>%
  filter(complete.cases(.))

#===============================================================================
# 3. SAVE CLEANED DATA
#===============================================================================

saveRDS(ukhls_complete_cases, file = paste0(OUTPUT_DIR, "ukhls_complete_cases.rds"))
