################################################################################
# Dissertation Project: The Effect of Smokers transitioning to E-cigarette in
#                       Physical and Mental health : An Emulated Trial using Longitudinal Data
#
# SCRIPT 3: PROPENSITY SCORE MATCHING (PSM)
#
# What this script does:
# 1. Loads clean data and creates the paired cohort.
# 2. Performs 1:3 (main) and 1:1 (sensitivity) propensity score matching.
# 3. Prints matching summaries to show sample size differences.
# 4. Assesses covariate balance for the main matching and generates Love plots.
# 5. Creates a Venn diagram to show overlap in matched individuals.
# 6. Saves the three final 1:3 matched datasets for the main analysis.
#
################################################################################


#===============================================================================
# 1. SETUP
#===============================================================================
library(tidyverse)
library(MatchIt)  # For propensity score matching
library(cobalt)   # For balance assessment and Love plots
library(ggvenn)   # For Venn diagrams

# --- Define paths ---
INPUT_DIR <- "./output_data/"
OUTPUT_DIR <- "./output_results/"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

#===============================================================================
# 2. LOAD DATA AND CREATE PAIRED COHORT
#===============================================================================

# --- Load the clean dataset from Script 1 ---
ukhls_paired <- readRDS(paste0(INPUT_DIR, "ukhls_paired.rds"))

#===============================================================================
# 3. PREPARE DATA AND DEFINITIONS FOR MATCHING
#===============================================================================

# --- Create three separate datasets for each comparison vs. 'Continued Smoker' ---
continue_vs_switcher <- filter(ukhls_paired, smoking_group %in% c("Continued Smoker", "Switcher")) %>%
  mutate(treatment = ifelse(smoking_group == "Switcher", 1, 0))

continue_vs_dual <- filter(ukhls_paired, smoking_group %in% c("Continued Smoker", "Dual User")) %>%
  mutate(treatment = ifelse(smoking_group == "Dual User", 1, 0))

continue_vs_quitter <- filter(ukhls_paired, smoking_group %in% c("Continued Smoker", "Quitter")) %>%
  mutate(treatment = ifelse(smoking_group == "Quitter", 1, 0))

# --- Define the propensity score model formula ---
match_formula <- as.formula(treatment ~ age_dv_t0 + age_squared_t0 + sex_dv_t0 + ethn_dv_t0 + hiqual_dv_t0 +
                              health_t0 + sf12pcs_dv_t0 + sf12mcs_dv_t0 + jbstat_t0 + log_real_income_t0 +
                              hl2gp_t0 + gor_dv_t0 + nkids015_t0 + factor(wave_t0))
# --- Define pretty names for variables for use in plots ---
pretty_variable_names <- c(
  age_dv_t0 = "Age", age_squared_t0 = "Age²", sex_dv_t0 = "Sex", ethn_dv_t0 = "Ethnicity",
  hiqual_dv_t0 = "Education", health_t0 = "Long-term Illness", jbstat_t0 = "Employment Status",
  log_real_income_t0 = "Log Real Income", sf12pcs_dv_t0 = "SF-12 PCS (Baseline)",
  sf12mcs_dv_t0 = "SF-12 MCS (Baseline)", hl2gp_t0 = "GP Visits", gor_dv_t0 = "Region",
  nkids015_t0 = "No. of Children", `factor(wave_t0)` = "Wave", distance = "Propensity Score"
)

#===============================================================================
# 4. PERFORM PROPENSITY SCORE MATCHING
#===============================================================================

# --- Main Analysis: 1:3 Nearest Neighbor Matching ---
match_switcher_1to3 <- matchit(match_formula, data = continue_vs_switcher, method = "nearest", ratio = 3)
match_dual_1to3     <- matchit(match_formula, data = continue_vs_dual,     method = "nearest", ratio = 3)
match_quitter_1to3  <- matchit(match_formula, data = continue_vs_quitter,  method = "nearest", ratio = 3)


# --- Sensitivity Analysis: 1:1 Nearest Neighbor Matching ---
match_switcher_1to1 <- matchit(match_formula, data = continue_vs_switcher, method = "nearest", ratio = 1)
match_dual_1to1     <- matchit(match_formula, data = continue_vs_dual,     method = "nearest", ratio = 1)
match_quitter_1to1  <- matchit(match_formula, data = continue_vs_quitter,  method = "nearest", ratio = 1)

#===============================================================================
# 5. ASSESS BALANCE AND VISUALIZE
#===============================================================================

# --- Generate and save Love plots for the main 1:3 matching analysis ---
love_plot_switcher <- love.plot(match_switcher_1to3, stat = "mean.diffs", abs = TRUE,
                                thresholds = c(m = 0.1), stars = "std", var.order = "unadjusted",
                                title = "Covariate Balance: Continued Smokers vs. Switchers (1:3 Matching)",
                                var.names = pretty_variable_names, wrap = 40) +
  theme(plot.title = element_text(hjust = 0.5))

love_plot_dual <- love.plot(match_dual_1to3, stat = "mean.diffs", abs = TRUE,
                            thresholds = c(m = 0.1), stars = "std", var.order = "unadjusted",
                            title = "Covariate Balance: Continued Smokers vs. Dual Users (1:3 Matching)",
                            var.names = pretty_variable_names, wrap = 40) +
  theme(plot.title = element_text(hjust = 0.5))

love_plot_quitter <- love.plot(match_quitter_1to3, stat = "mean.diffs", abs = TRUE,
                               thresholds = c(m = 0.1), stars = "std", var.order = "unadjusted",
                               title = "Covariate Balance: Continued Smokers vs. Quitters (1:3 Matching)",
                               var.names = pretty_variable_names, wrap = 40) +
  theme(plot.title = element_text(hjust = 0.5))

# Save plots
ggsave(paste0(OUTPUT_DIR, "LovePlot_vs_Switcher_1to3.png"), plot = love_plot_switcher, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_DIR, "LovePlot_vs_DualUser_1to3.png"), plot = love_plot_dual, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(paste0(OUTPUT_DIR, "LovePlot_vs_Quitter_1to3.png"), plot = love_plot_quitter, width = 10, height = 7, dpi = 300, bg = "white")

# --- Extract and save the 1:3 matched dataframes for subsequent analyses ---
matched_data_switcher <- match.data(match_switcher_1to3)
matched_data_dual     <- match.data(match_dual_1to3)
matched_data_quitter  <- match.data(match_quitter_1to3)

# --- Create Venn Diagram of overlapped individuals from the main 1:3 matched data ---
venn_list <- list(
  `vs. Switcher` = matched_data_switcher %>%
    mutate(unique_id = paste(pidp, wave_pair, sep = "_")) %>%
    pull(unique_id),
  
  `vs. Dual User` = matched_data_dual %>%
    mutate(unique_id = paste(pidp, wave_pair, sep = "_")) %>%
    pull(unique_id),
  
  `vs. Quitter` = matched_data_quitter %>%
    mutate(unique_id = paste(pidp, wave_pair, sep = "_")) %>%
    pull(unique_id)
)

venn_diagram <- ggvenn(
  venn_list,
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  set_name_size = 5, text_size = 4
) +
  labs(
    title = "Overlap of Unique Individuals Across 1:3 Matched Analyses",
    caption = "Each circle represents a full matched dataset (Controls + Treatment)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

print(venn_diagram)
ggsave(paste0(OUTPUT_DIR, "VennDiagram_Matched_Overlap_1to3.png"), plot = venn_diagram, width = 8, height = 8, dpi = 300, bg = "white")


#===============================================================================
# 6. SAVE FINAL MATCHED DATASETS
#===============================================================================

saveRDS(matched_data_switcher, file = paste0(INPUT_DIR, "matched_data_switcher_1to3.rds"))
saveRDS(matched_data_dual,     file = paste0(INPUT_DIR, "matched_data_dual_1to3.rds"))
saveRDS(matched_data_quitter,  file = paste0(INPUT_DIR, "matched_data_quitter_1to3.rds"))
