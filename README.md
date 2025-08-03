The Effect of Smokers Transitioning to E-cigarettes on Physical and Mental Health: An Emulated Trial using Longitudinal Data

Project Overview

This repository contains the full R code and analysis for a dissertation project investigating the health impacts of smokers switching to e-cigarettes. The project emulates a target trial using observational data from the UK Household Longitudinal Study (UKHLS) to estimate the causal effect of this transition on physical health (SF-12 PCS), mental health (SF-12 MCS), and overall health-related quality of life (EQ-5D).

The analysis follows a rigorous, multi-stage approach, including:

Propensity Score Matching (PSM) to create comparable groups of "Switchers" and "Continued Smokers".

Doubly Robust Regression Models to control for confounding variables.

Fixed-Effects Models as a sensitivity analysis to account for time-invariant individual characteristics.

Subgroup Analyses to explore effects across different demographic groups.
Of course. A good README file is essential for any project. It acts as the front page, explaining what the project is, why it's important, and how to use it.

Here is a comprehensive README file for your dissertation project, written in markdown format. You can copy and paste this directly into a new file named README.md in the main directory of your GitHub repository.

The Effect of Smokers Transitioning to E-cigarettes on Physical and Mental Health
An Emulated Trial using Longitudinal Data from the UK Household Longitudinal Study (UKHLS)
📖 Project Overview
This repository contains the full R code and analysis for a dissertation project investigating the health impacts of smokers switching to e-cigarettes. The project emulates a target trial using observational data from the UK Household Longitudinal Study (UKHLS) to estimate the causal effect of this transition on physical health (SF-12 PCS), mental health (SF-12 MCS), and overall health-related quality of life (EQ-5D).

The analysis follows a rigorous, multi-stage approach, including:

Propensity Score Matching (PSM) to create comparable groups of "Switchers" and "Continued Smokers".

Doubly Robust Regression Models to control for confounding variables.

Fixed-Effects Models as a sensitivity analysis to account for time-invariant individual characteristics.

Subgroup Analyses to explore effects across different demographic groups.

Project Structure
The project is organized into a sequential pipeline of R scripts, designed to be run in order. The main directories are:

/output_data/: Stores intermediate, cleaned datasets (.rds files) that are passed between scripts.

/output_results/: Stores all numerical results (.rds files) and diagnostic plots.

/output_results/final_tables_and_figures/: Stores all final, publication-quality tables and figures (.png, .docx).

Scripts Pipeline
The analysis is broken down into a series of modular scripts:

1_data_preparation.R: Loads the raw Stata dataset, performs all cleaning and recoding, creates the primary "paired" analytical dataset, and saves the clean data to the /output_data/ folder.

2_descriptive_and_regression.R: Loads the clean paired data and conducts the initial analysis on the full unmatched cohort. It generates a baseline characteristics table and runs unadjusted/adjusted regression models.

3_propensity_score_matching.R: Performs the core propensity score matching (both 1:1 and 1:3 ratios). It assesses covariate balance by generating tables and Love plots.

4_matched_regression_analysis.R: Runs the main analysis on the matched cohort. This includes both non-doubly robust and doubly robust regression models.

5_sensitivity_and_subgroup.R: Conducts the final analytical steps, including Fixed-Effects models as a sensitivity analysis and subgroup analyses by age, education, and income.

6_final_reporting.R: Loads all previously saved results and generates every final, publication-quality table and figure for the dissertation appendix.
