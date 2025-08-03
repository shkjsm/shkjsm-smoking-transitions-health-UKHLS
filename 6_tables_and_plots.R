################################################################################
# Dissertation Project: The Effect of Smokers transitioning to E-cigarette in
#                       Physical and Mental health : An Emulated Trial using Longitudinal Data
#
# SCRIPT 5: SENSITIVITY AND SUBGROUP ANALYSIS
#
# What this script does:
# 1. Loads the 1:3 matched datasets.
# 2. Runs Fixed-Effects (FE) models as a sensitivity analysis.
# 3. Runs doubly robust models stratified by age, education, and income.
# 4. Generates and saves final summary tables for each analysis.
#
################################################################################


#===============================================================================
# 1. SETUP
#===============================================================================
library(tidyverse)
library(plm)      # For Fixed-Effects models
library(lmtest)   # For robust standard errors (coeftest)
library(estimatr) # For lm_robust in subgroup models
library(broom)
library(gt)
library(officer) # To save table as .docx

# --- Define paths ---
INPUT_DIR <- "./output_data/"
OUTPUT_DIR <- "./output_results/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

#===============================================================================
# 2. LOAD AND PREPARE DATA
#===============================================================================

# --- Load the 1:3 matched datasets ---
matched_switcher <- readRDS(paste0(INPUT_DIR, "matched_data_switcher_1to3.rds"))
matched_dual     <- readRDS(paste0(INPUT_DIR, "matched_data_dual_1to3.rds"))
matched_quitter  <- readRDS(paste0(INPUT_DIR, "matched_data_quitter_1to3.rds"))

# --- Helper function to prepare data (creates subgroup vars) ---
prepare_analysis_data <- function(data) {
  data %>%
    mutate(
      age_group_t0 = case_when(
        age_dv_t0 >= 16 & age_dv_t0 <= 24 ~ "16-24",
        age_dv_t0 >= 25 & age_dv_t0 <= 49 ~ "25-49",
        age_dv_t0 >= 50 ~ "50+"),
      income_quintile_t0 = factor(ntile(log_real_income_t0, 5))
    )
}

data_switcher <- prepare_analysis_data(matched_switcher)
data_dual     <- prepare_analysis_data(matched_dual)
data_quitter  <- prepare_analysis_data(matched_quitter)

#===============================================================================
# 3. SENSITIVITY ANALYSIS: FIXED-EFFECTS MODELS
#===============================================================================

# --- Define FE Model Formulas (time-varying covariates only) ---
fe_formula_pcs <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + age_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + nkids015_t0)
fe_formula_mcs <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12mcs_dv_t0 + age_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + nkids015_t0)
fe_formula_eq5d <- as.formula(eq5d_t1 ~ smoking_group + eq5d_t0 + age_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + nkids015_t0)

# --- Reusable function to run an FE model ---
run_fe_model <- function(formula, data) {
  pdata <- pdata.frame(data, index = c("pidp", "wave_t1"))
  model <- plm(formula, data = pdata, model = "within")
  # Get robust standard errors clustered by individual
  coeftest(model, vcov = vcovHC(model, type = "HC1", cluster = "group")) %>%
    tidy(conf.int = TRUE)
}

# --- Run all 9 FE Models ---

FE_PCS_vs_Switcher = run_fe_model(fe_formula_pcs, data_switcher)
FE_MCS_vs_Switcher = run_fe_model(fe_formula_mcs, data_switcher)
FE_eq5d_vs_Switcher = run_fe_model(fe_formula_eq5d, data_switcher)
FE_PCS_vs_Dual = run_fe_model(fe_formula_pcs, data_dual)
FE_MCS_vs_Dual = run_fe_model(fe_formula_mcs, data_dual)
FE_eq5d_vs_Dual = run_fe_model(fe_formula_eq5d, data_dual)
FE_PCS_vs_Quitter = run_fe_model(fe_formula_pcs, data_quitter)
FE_MCS_vs_Quitter = run_fe_model(fe_formula_mcs, data_quitter)
FE_eq5d_vs_Quitter = run_fe_model(fe_formula_eq5d, data_quitter)

#---------- save as list -----------------
fe_results_list <- list(
  FE_PCS_vs_Switcher = run_fe_model(fe_formula_pcs, data_switcher),
  FE_MCS_vs_Switcher = run_fe_model(fe_formula_mcs, data_switcher),
  FE_eq5d_vs_Switcher = run_fe_model(fe_formula_eq5d, data_switcher),
  FE_PCS_vs_Dual = run_fe_model(fe_formula_pcs, data_dual),
  FE_MCS_vs_Dual = run_fe_model(fe_formula_mcs, data_dual),
  FE_eq5d_vs_Dual = run_fe_model(fe_formula_eq5d, data_dual),
  FE_PCS_vs_Quitter = run_fe_model(fe_formula_pcs, data_quitter),
  FE_MCS_vs_Quitter = run_fe_model(fe_formula_mcs, data_quitter),
  FE_eq5d_vs_Quitter = run_fe_model(fe_formula_eq5d, data_quitter)
)

#===============================================================================
# 4. SUBGROUP ANALYSES
#===============================================================================

# --- Define DR model formula (same as Script 4) ---
pcs_formula <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + wave_t0)
mcs_formula <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + wave_t0)
pcs_formula_no_educ <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + wave_t0)
mcs_formula_no_educ <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + wave_t0)

# --- Reusable function to run subgroup models ---
run_subgroup_model <- function(data, group_var, formula) {
  data %>%
    group_by({{ group_var }}) %>%
    group_modify(~ tidy(lm_robust(formula, data = .x, clusters = pidp))) %>%
    ungroup()
}
# --- Run all subgroup models ---
subgroup_results_list <- list(
  Age_PCS_vs_Switcher = run_subgroup_model(data_switcher, age_group_t0, pcs_formula),
  Age_PCS_vs_Dual = run_subgroup_model(data_dual, age_group_t0, pcs_formula),
  Age_PCS_vs_Quitter = run_subgroup_model(data_quitter, age_group_t0, pcs_formula),
  Age_MCS_vs_Switcher = run_subgroup_model(data_switcher, age_group_t0, mcs_formula),
  Age_MCS_vs_Dual = run_subgroup_model(data_dual, age_group_t0, mcs_formula),
  Age_MCS_vs_Quitter = run_subgroup_model(data_quitter, age_group_t0, mcs_formula),
  Educ_PCS_vs_Switcher = run_subgroup_model(data_switcher, hiqual_dv_t0, pcs_formula_no_educ),
  Educ_PCS_vs_Dual = run_subgroup_model(data_dual, hiqual_dv_t0, pcs_formula_no_educ),
  Educ_PCS_vs_Quitter = run_subgroup_model(data_quitter, hiqual_dv_t0, pcs_formula_no_educ),
  Educ_MCS_vs_Switcher = run_subgroup_model(data_switcher, hiqual_dv_t0, mcs_formula_no_educ),
  Educ_MCS_vs_Dual = run_subgroup_model(data_dual, hiqual_dv_t0, mcs_formula_no_educ),
  Educ_MCS_vs_Quitter = run_subgroup_model(data_quitter, hiqual_dv_t0, mcs_formula_no_educ),
  Income_PCS_vs_Switcher = run_subgroup_model(data_switcher, income_quintile_t0, pcs_formula),
  Income_PCS_vs_Dual = run_subgroup_model(data_dual, income_quintile_t0, pcs_formula),
  Income_PCS_vs_Quitter = run_subgroup_model(data_quitter, income_quintile_t0, pcs_formula),
  Income_MCS_vs_Switcher = run_subgroup_model(data_switcher, income_quintile_t0, mcs_formula),
  Income_MCS_vs_Dual = run_subgroup_model(data_dual, income_quintile_t0, mcs_formula),
  Income_MCS_vs_Quitter = run_subgroup_model(data_quitter, income_quintile_t0, mcs_formula)
)

#===============================================================================
# 5. GENERATE AND SAVE FINAL TABLES
#===============================================================================

# --- Helper function for formatting ---
format_results <- function(estimate, conf.low, conf.high, p.value) {
  stars <- case_when(p.value < 0.001 ~ "***", p.value < 0.01  ~ "**", p.value < 0.05  ~ "*", TRUE ~ "")
  p_formatted <- case_when(p.value < 0.001 ~ "*p*<0.001", TRUE ~ paste0("*p*=", sprintf("%.3f", p.value)))
  paste0(sprintf("%.3f", estimate), stars, " [", sprintf("%.3f", conf.low), ", ", sprintf("%.3f", conf.high), "], ", p_formatted)
}

# --- Prepare Data for the Table ---
fe_results_data <- bind_rows(
  FE_PCS_vs_Switcher %>% mutate(outcome = "Physical Health (SF-12 PCS)"),
  FE_MCS_vs_Switcher %>% mutate(outcome = "Mental Health (SF-12 MCS)"),
  FE_eq5d_vs_Switcher %>% mutate(outcome = "Health Related Quality of Life (EQ-5D-3L)"),
  FE_PCS_vs_Dual %>% mutate(outcome = "Physical Health (SF-12 PCS)"),
  FE_MCS_vs_Dual %>% mutate(outcome = "Mental Health (SF-12 MCS)"),
  FE_eq5d_vs_Dual %>% mutate(outcome = "Health Related Quality of Life (EQ-5D-3L)"),
  FE_PCS_vs_Quitter %>% mutate(outcome = "Physical Health (SF-12 PCS)"),
  FE_MCS_vs_Quitter %>% mutate(outcome = "Mental Health (SF-12 MCS)"),
  FE_eq5d_vs_Quitter%>% mutate(outcome = "Health Related Quality of Life (EQ-5D-3L)")
) %>%
  filter(str_detect(term, "smoking_group")) %>%
  mutate(
    Comparison = str_remove(term, "smoking_group"),
    Result = format_results(estimate, conf.low, conf.high, p.value)
  ) %>%
  select(outcome, Comparison, Result)

# --- Create the Final Table with gt ---
fe_results_gt <- gt(fe_table_data, groupname_col = "Outcome", rowname_col = "Comparison") %>%
  # Add titles
  tab_header(
    title = md("**Sensitivity Analysis of Smoking Transitions on Health outcomes**"),
    subtitle = "Results from Fixed-Effects Models on the Matched Cohort"
  ) %>%
  # Set the single column label
  cols_label(
    Result = md("**Est. [95% CI], *p* value**")
  ) %>%
  # Center align all columns
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  # Style the row group labels to be bold
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  # Render the markdown for italics and bold stars
  fmt_markdown(columns = everything()) %>%
  # Add footnotes
  tab_footnote(
    footnote = md("*Note: Reference group is 'Continued Smoker'*. <br>**p<0.05*, ***p<0.01*, ****p<0.001*")
  ) %>%
  # Final styling
  tab_options(table.width = pct(100))

# Print the table to the viewer
print(fe_results_gt)

# Save the table as a high-quality PNG
gtsave(fe_results_gt, filename = paste0(OUTPUT_DIR, "Table_Fixed_Effects_Results_Formatted.png"))

#save to generate comprehensive results later
saveRDS(fe_results_data, file = paste0(INPUT_DIR, "FE_Results"))

# --- Table 5: Subgroup Analysis Results ---
subgroup_table_data <- subgroup_results_list %>%
  bind_rows(.id="id") %>%
  filter(str_detect(term, "smoking_group")) %>%
  separate(id, into=c("Subgroup", "outcome", "temp", "Comparison"), sep="_") %>%
  select(-temp) %>%
  mutate(
    Result = format_results(estimate, conf.low, conf.high, p.value),
    Subgroup_Value = coalesce(age_group_t0, hiqual_dv_t0, as.character(income_quintile_t0))
  ) %>%
  select(Subgroup, Subgroup_Value, outcome, Comparison, Result) %>%
  pivot_wider(names_from = c(outcome, Comparison), values_from = Result)

subgroup_gt <- gt(subgroup_table_data, groupname_col = "Subgroup", rowname_col = "Subgroup_Value") %>%
  tab_header(title = md("**Subgroup Analysis for Health outcomes**"),
             subtitle = "Estimates from Doubly Robust Models") %>%
  fmt_markdown(columns = everything()) %>%
  tab_footnote(footnote = md("Note: Reference group is 'Continued Smoker'. *p<0.05, **p<0.01, ***p<0.001"))

gtsave(subgroup_gt, filename = paste0(OUTPUT_DIR, "Table5_Subgroup_Analysis_Results.docx"))

#save to geenrate comprehensive results later
saveRDS(subgroup_gt, file = paste0(INPUT_DIR, "Subgroup_Results"))
