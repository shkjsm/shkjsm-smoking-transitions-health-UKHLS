################################################################################
# Dissertation Project: The Effect of Smokers transitioning to E-cigarette in
#                       Physical and Mental health : An Emulated Trial using Longitudinal Data
#
# SCRIPT 4: REGRESSION ANALYSIS ON MATCHED DATA
#
# What this script does:
# 1. Loads the three 1:3 matched datasets.
# 2. Maps EQ-5D-3L scores for the matched data.
# 3. Runs non-doubly robust and doubly robust models.
# 4. Performs and saves detailed model diagnostics for the main models.
# 5. Generates and saves a final summary table as a PNG image.
#
################################################################################


#===============================================================================
# 1. SETUP
#===============================================================================
library(tidyverse)
library(estimatr) # For robust standard errors
library(broom)    # For tidying model output
library(gt)       # For publication-quality tables
library(car)      # For VIF (diagnostics)
library(lmtest)   # For Breusch-Pagan test (diagnostics)

# --- Define paths ---
INPUT_DIR <- "./output_data/"
OUTPUT_DIR <- "./output_results/"

#===============================================================================
# 2. LOAD AND PREPARE MATCHED DATA
#===============================================================================

# --- Load the 1:3 matched datasets from Script 3 ---
matched_data_switcher <- readRDS(paste0(INPUT_DIR, "matched_data_switcher_1to3.rds"))
matched_data_dual     <- readRDS(paste0(INPUT_DIR, "matched_data_dual_1to3.rds"))
matched_data_quitter  <- readRDS(paste0(INPUT_DIR, "matched_data_quitter_1to3.rds"))

#===============================================================================
# 3. REGRESSION MODEL ESTIMATION
#===============================================================================

# --- Define Model Equations ---
pcs_formula <- as.formula(sf12pcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
mcs_formula <- as.formula(sf12mcs_dv_t1 ~ smoking_group + sf12pcs_dv_t0 + sf12mcs_dv_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
eq5d_formula <- as.formula(eq5d_t1 ~ smoking_group + eq5d_t0 + age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 + health_t0 + jbstat_t0 + log_real_income_t0 + hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))

# --- Non-Doubly Robust Models (Unadjusted on Matched Data) ---
ndr_switcher_pcs  <- tidy(lm_robust(sf12pcs_dv_t1 ~ smoking_group, data = matched_data_switcher, clusters = pidp))
ndr_switcher_mcs  <- tidy(lm_robust(sf12mcs_dv_t1 ~ smoking_group, data = matched_data_switcher, clusters = pidp))
ndr_switcher_eq5d <- tidy(lm_robust(eq5d_t1 ~ smoking_group,         data = matched_data_switcher, clusters = pidp))

ndr_dual_pcs      <- tidy(lm_robust(sf12pcs_dv_t1 ~ smoking_group, data = matched_data_dual, clusters = pidp))
ndr_dual_mcs      <- tidy(lm_robust(sf12mcs_dv_t1 ~ smoking_group, data = matched_data_dual, clusters = pidp))
ndr_dual_eq5d     <- tidy(lm_robust(eq5d_t1 ~ smoking_group,         data = matched_data_dual, clusters = pidp))

ndr_quitter_pcs   <- tidy(lm_robust(sf12pcs_dv_t1 ~ smoking_group, data = matched_data_quitter, clusters = pidp))
ndr_quitter_mcs   <- tidy(lm_robust(sf12mcs_dv_t1 ~ smoking_group, data = matched_data_quitter, clusters = pidp))
ndr_quitter_eq5d  <- tidy(lm_robust(eq5d_t1 ~ smoking_group,         data = matched_data_quitter, clusters = pidp))

# --- Doubly Robust Models (Adjusted on Matched Data) ---
dr_switcher_pcs  <- tidy(lm_robust(pcs_formula, data = matched_data_switcher, clusters = pidp))
dr_switcher_mcs  <- tidy(lm_robust(mcs_formula, data = matched_data_switcher, clusters = pidp))
dr_switcher_eq5d <- tidy(lm_robust(eq5d_formula, data = matched_data_switcher, clusters = pidp))

dr_dual_pcs      <- tidy(lm_robust(pcs_formula, data = matched_data_dual, clusters = pidp))
dr_dual_mcs      <- tidy(lm_robust(mcs_formula, data = matched_data_dual, clusters = pidp))
dr_dual_eq5d     <- tidy(lm_robust(eq5d_formula, data = matched_data_dual, clusters = pidp))

dr_quitter_pcs   <- tidy(lm_robust(pcs_formula, data = matched_data_quitter, clusters = pidp))
dr_quitter_mcs   <- tidy(lm_robust(mcs_formula, data = matched_data_quitter, clusters = pidp))
dr_quitter_eq5d  <- tidy(lm_robust(eq5d_formula, data = matched_data_quitter, clusters = pidp))

#===============================================================================
# 4. MODEL DIAGNOSTICS FOR DOUBLY ROBUST MODELS
#===============================================================================

# --- Re-run DR models with standard lm() to access diagnostic functions ---
DR_Switcher_PCS = lm(pcs_formula, data = matched_data_switcher)
DR_Switcher_MCS = lm(mcs_formula, data = matched_data_switcher)
DR_Switcher_eq5d = lm(eq5d_formula, data = matched_data_switcher)
DR_Dual_PCS = lm(pcs_formula, data = matched_data_dual)
DR_Dual_MCS = lm(mcs_formula, data = matched_data_dual)
DR_Dual_eq5d = lm(eq5d_formula, data = matched_data_dual)
DR_Quitter_PCS = lm(pcs_formula, data = matched_data_quitter)
DR_Quitter_MCS = lm(mcs_formula, data = matched_data_quitter)
DR_Quitter_eq5d = lm(eq5d_formula, data = matched_data_quitter)

# --- Plots for SF-12 PCS Models ---
png("Diagnostic_Plots_PCS_vs_Switcher.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Switcher_PCS, main = "Diagnostics: PCS vs. Switcher")
dev.off()

png("Diagnostic_Plots_PCS_vs_DualUser.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Dual_PCS, main = "Diagnostics: PCS vs. Dual User")
dev.off()

png("Diagnostic_Plots_PCS_vs_Quitter.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Quitter_PCS, main = "Diagnostics: PCS vs. Quitter")
dev.off()


# --- Plots for SF-12 MCS Models ---
png("Diagnostic_Plots_MCS_vs_Switcher.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Switcher_MCS, main = "Diagnostics: MCS vs. Switcher")
dev.off()

png("Diagnostic_Plots_MCS_vs_DualUser.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Dual_MCS, main = "Diagnostics: MCS vs. Dual User")
dev.off()

png("Diagnostic_Plots_MCS_vs_Quitter.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Quitter_MCS, main = "Diagnostics: MCS vs. Quitter")
dev.off()


# --- Plots for EQ-5D-3L Models ---
png("Diagnostic_Plots_EQ5D_vs_Switcher.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Switcher_eq5d, main = "Diagnostics: EQ-5D-3L vs. Switcher")
dev.off()

png("Diagnostic_Plots_EQ5D_vs_DualUser.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Dual_eq5d, main = "Diagnostics: EQ-5D-3L vs. Dual User")
dev.off()

png("Diagnostic_Plots_EQ5D_vs_Quitter.png", width = 1000, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(DR_Quitter_eq5d, main = "Diagnostics: EQ-5D-3L vs. Quitter")
dev.off()

#===============================================================================
# 5. GENERATE AND SAVE FINAL RESULTS TABLE
#===============================================================================

# --- Helper function for formatting result strings ---
format_results <- function(estimate, conf.low, conf.high, p.value) {
  stars <- case_when(p.value < 0.001 ~ "***", p.value < 0.01  ~ "**", p.value < 0.05  ~ "*", TRUE ~ "")
  p_formatted <- case_when(p.value < 0.001 ~ "*p*<0.001", TRUE ~ paste0("*p*=", sprintf("%.3f", p.value)))
  paste0(sprintf("%.3f", estimate), stars, " [", sprintf("%.3f", conf.low), ", ", sprintf("%.3f", conf.high), "], ", p_formatted)
}

# --- Combine all results into a single dataframe ---
matched_results_data <- bind_rows(
  "Physical Health (SF-12 PCS)" = bind_rows(
    "Unadjusted (Non-Doubly Robust)" = bind_rows("Switcher" = ndr_switcher_pcs, "Dual User" = ndr_dual_pcs, "Quitter" = ndr_quitter_pcs, .id="Comparison"),
    "Fully Adjusted (Doubly Robust)" = bind_rows("Switcher" = dr_switcher_pcs, "Dual User" = dr_dual_pcs, "Quitter" = dr_quitter_pcs, .id="Comparison"), .id="Model"),
  "Mental Health (SF-12 MCS)" = bind_rows(
    "Unadjusted (Non-Doubly Robust)" = bind_rows("Switcher" = ndr_switcher_mcs, "Dual User" = ndr_dual_mcs, "Quitter" = ndr_quitter_mcs, .id="Comparison"),
    "Fully Adjusted (Doubly Robust)" = bind_rows("Switcher" = dr_switcher_mcs, "Dual User" = dr_dual_mcs, "Quitter" = dr_quitter_mcs, .id="Comparison"), .id="Model"),
  "Health Related Quality of Life (EQ-5D-3L)" = bind_rows(
    "Unadjusted (Non-Doubly Robust)" = bind_rows("Switcher" = ndr_switcher_eq5d, "Dual User" = ndr_dual_eq5d, "Quitter" = ndr_quitter_eq5d, .id="Comparison"),
    "Fully Adjusted (Doubly Robust)" = bind_rows("Switcher" = dr_switcher_eq5d, "Dual User" = dr_dual_eq5d, "Quitter" = dr_quitter_eq5d, .id="Comparison"), .id="Model"),
  .id = "outcome"
) %>%
  filter(str_detect(term, "smoking_group")) %>%
  mutate(Result = format_results(estimate, conf.low, conf.high, p.value)) %>%
  select(outcome, Comparison, Model, Result) %>%
  pivot_wider(names_from = Model, values_from = Result)

# --- Create and save the final table ---
matched_results_gt <- gt(matched_results_data, groupname_col = "outcome", rowname_col = "Comparison") %>%
  tab_header(
    title = md("**Associations Between Smoking Transitions and Health outcomes (1:3 PS matched Cohort)**"),
    subtitle = "Results from Models on the Matched Cohort"
  ) %>%
  cols_label(
    `Unadjusted (Non-Doubly Robust)` = md("**Unadjusted (Matched)**<br>Est. [95% CI], *p*-value"),
    `Fully Adjusted (Doubly Robust)` = md("**Fully Adjusted (Dobuly Robust)**<br>Est. [95% CI], *p*-value")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups())%>%
  fmt_markdown(columns = everything()) %>%
  tab_options(table.width = pct(100)) %>%
  tab_footnote(footnote = md("*Note: Reference group is 'Continued Smoker'*.<br>**p*<0.05, ***p*<0.01, ****p*<0.001"))

print(matched_results_gt)
gtsave(matched_results_gt, filename = paste0(OUTPUT_DIR, "Table3_Matched_Results.png"))

#save to generate comprehensice results later
saveRDS(matched_results_data, file = paste0(INPUT_DIR, "Matched_Results"))


##############---------------------------------------------#####################
#------------- Full regression Table Generation --------------------------------
##############---------------------------------------------#####################

# --- Add model metadata before binding ---
full_output <- function(model_df, group, outcome) {
  model_df %>%
    mutate(
      SmokingGroup = group,
      Outcome = outcome
    )
}

# Combine all DR models with metadata
full_dr_results <- bind_rows(
  full_output(dr_switcher_pcs, "Switcher", "SF-12 PCS"),
  full_output(dr_dual_pcs,     "Dual User", "SF-12 PCS"),
  full_output(dr_quitter_pcs,  "Quitter",   "SF-12 PCS"),
  
  full_output(dr_switcher_mcs, "Switcher", "SF-12 MCS"),
  full_output(dr_dual_mcs,     "Dual User", "SF-12 MCS"),
  full_output(dr_quitter_mcs,  "Quitter",   "SF-12 MCS"),
  
  full_output(dr_switcher_eq5d,"Switcher", "EQ-5D-3L"),
  full_output(dr_dual_eq5d,    "Dual User", "EQ-5D-3L"),
  full_output(dr_quitter_eq5d, "Quitter",   "EQ-5D-3L")
)

# Optional: format term labels
clean_term_names <- function(term) {
  term %>%
    str_replace_all("smoking_group", "Smoking Group: ") %>%
    str_replace_all("sf12pcs_dv_t0", "Baseline PCS") %>%
    str_replace_all("sf12mcs_dv_t0", "Baseline MCS") %>%
    str_replace_all("eq5d_t0", "Baseline EQ-5D-3L") %>%
    str_replace_all("age_dv_t0", "Age") %>%
    str_replace_all("age_squared_t0", "Age²") %>%
    str_replace_all("sex_dv_t0", "Sex: ") %>%
    str_replace_all("ethn_dv_t0", "Ethnicity: ") %>%
    str_replace_all("hiqual_dv_t0", "Highest Education: ") %>%
    str_replace_all("health_t0", "Longterm Illness: ") %>%
    str_replace_all("jbstat_t0", "Employment Status: ") %>%
    str_replace_all("log_real_income_t0", "Log Income: ") %>%
    str_replace_all("hl2gp_t0", "Annual GP Visits: ") %>%
    str_replace_all("gor_dv_t0", "Region: ") %>%
    str_replace_all("nkids015_t0", "No. of Kids: ") %>%
    str_replace_all("wave_t0", "Survey Wave")%>%
    str_replace_all("sf12pcs_dv_t1", "Physical Health (SF-12 PCS)")%>%
    str_replace_all("sf12mcs_dv_t1", "Mental Health (SF-12 MCS)")%>%
    str_replace_all("eq-5d_t1", "Health Related Quality of Life (EQ-5D-3L")
}

# Add formatted result column
format_results <- function(estimate, conf.low, conf.high, p.value) {
  stars <- case_when(p.value < 0.001 ~ "***", p.value < 0.01  ~ "**", p.value < 0.05  ~ "*", TRUE ~ "")
  p_formatted <- case_when(p.value < 0.001 ~ "*p*<0.001", TRUE ~ paste0("*p*=", sprintf("%.3f", p.value)))
  paste0(sprintf("%.3f", estimate), stars, " [", sprintf("%.3f", conf.low), ", ", sprintf("%.3f", conf.high), "], ", p_formatted)
}
full_dr_results_formatted <- full_dr_results %>%
  mutate(
    term = clean_term_names(term),
    Result = format_results(estimate, conf.low, conf.high, p.value)
  ) %>%
  select(Outcome, SmokingGroup, term, Result)

# Create wide table: terms as rows, outcomes and groups as columns
regression_table_wide <- full_dr_results_formatted %>%
  pivot_wider(names_from = SmokingGroup, values_from = Result)

#create table
regression_table_gt <- regression_table_wide %>%
  gt(groupname_col = "Outcome", rowname_col = "term") %>%
  tab_header(
    title = md("**Full Regression Coefficients from Doubly Robust Models**"),
    subtitle = md("Outcomes: SF-12 PCS, SF-12 MCS, and EQ-5D")
  ) %>%
  cols_label(
    Switcher = md("**Switcher**<br>Est. [95% CI], *p*"),
    `Dual User` = md("**Dual User**<br>Est. [95% CI], *p*"),
    Quitter = md("**Quitter**<br>Est. [95% CI], *p*")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_options(table.width = pct(100))

print(regression_table_gt)

# Save the table
gtsave(regression_table_gt, filename = paste0(OUTPUT_DIR, "Full_DR_Regression_Table.docx"))
saveRDS(regression_table_wide, file = paste0(INPUT_DIR, "Full_DR_Regression_Table.rds"))


