# Correcting-for-effect-modification-in-the-doubly-ranked-non-linear-Mendelian-randomization-method

This directory contains all the scripts used in the manuscript on **correcting for effect modification in the doubly-ranked non-linear Mendelian randomization method**.

## Simulation Scenario 1

- **Description:** Scenario 1 examines how the magnitude of the net GxE interaction affects bias in the doubly-ranked method. The net interaction is varied by changing the correlation between two effect modifiers with opposing GxE interaction effects.
- **Data generation script:** `01_dgm_scenario1_20251009.R`
- **Data analysis script:** `02_analysis_scenario1_20251009.R`
- **Figure:** `Figure2.R`

## Simulation Scenario 2

- **Description:** Scenario 2 evaluates whether the proposed GxE correction can attenuate bias under null, linear, and non-linear exposure-outcome relationships.
- **Data generation script:** `03_dgm_scenario2_20251016.R`
- **Data analysis script:** `04_analysis_scenario2_20251016.R`
- **Figure:**
  - `Figure3.R`
  - `FigureS2.R` (mean-squared error for the quadratic setting)
  - `FigureS3.R` (results for the additional scenario with a threshold-effect of exposure on outcome)

## Simulation Scenario 3

- **Description:** Scenario 3 tests whether the GxE correction itself can introduce bias when the effect modifier also acts as a confounder or is highly correlated with a mediator or collider of the exposure-outcome association or itself acts as a collider.
- **Data generation script:**
  - `05_dgm_scenario3_20251016.R`
  - `07_dgm_scenarioSupp_20260330.R` (effect modifier is also a collider)
- **Data analysis script:**
  - `06_analysis_scenario3_20251016.R`
  - `08_analysis_scenarioSupp_20260330.R` (effect modifier is also a collider)
- **Figure:**
  - `Figure4.R`
  - `FigureS6.R` (effect modifier is also a collider)

## Supplementary Simulation Scenario 1

- **Description:** Supplementary scenario 1 evaluates the robustness of the correction when the true effect modifier is measured imperfectly.
- **Data generation script:** `07_dgm_scenarioSupp_20260330.R`
- **Data analysis script:** `08_analysis_scenarioSupp_20260330.R`
- **Figure:** `FigureS4.R`

## Supplementary Simulation Scenario 2

- **Description:** Supplementary scenario 2 evaluates the correction when the modifier is transformed onto a different scale.
- **Data generation script:** `07_dgm_scenarioSupp_20260330.R`
- **Data analysis script:** `08_analysis_scenarioSupp_20260330.R`
- **Figure:** `FigureS5.R`

## GxE Interaction Analyses

- **Description:** These scripts assess effect modification of the genetic score-exposure association across selected modifiers.
- **Data analysis scripts:**
  - `09_gxeInteraction_vitd_replicatedScore_20250405.R`
  - `10_gxeInteraction_vitd_focusedScore_20250405.R`
  - `11_gxeInteraction_bmi_20250407.R`
  - `12_gxeInteraction_ldl_20250428.R`
- **Figure:**
  - `Figure5.R`

## Falsification Test

- **Description:** Falsification test for vitamin D (using replicated score and focused score), BMI, and LDL-C.
- **Helper script:** `helpers/sim_nulloutcome.R`
- **Data analysis scripts:**
  - `13_falsification_vitd_replicatedScore_20251105.R`
  - `14_falsification_vitd_focusedScore_20251105.R`
  - `15_falsification_bmi_20251206.R`
  - `16_falsification_ldl_20251105.R`
- **Figure:**
  - `Figure6.R` (top two GxE corrections)
  - `Figure7.R` (normalized mean-squared error)
  - `FigureS7.R` (full results)

## Empirical Negative and Positive Control Analysis for LDL-C

- **Description:** This analysis follows up the falsification results for LDL-C using age and sex as empirical negative controls and CAD as a positive control.
- **Data analysis script:** `19_ldl_EmpNegPosControls_20251105.R`
- **Figure:**
  - `Figure8.R`
  - `FigureS10.R` (normalized mean-squared error across negative and positive control outcomes)

## Empirical Negative Control Analysis for Vitamin D Using Replicated Score

- **Description:** Empirical negative control analysis with age and sex for vitamin D using replicated score.
- **Data analysis script:** `17_vitd_EmpNegControls_replicatedScore_20251205.R`
- **Figure:** `FigureS8.R`

## Empirical Negative Control Analysis for Vitamin D Using Focused Score

- **Description:** Empirical negative control analysis with age and sex for vitamin D using focused score.
- **Data analysis:** `18_vitd_EmpNegControls_focusedScore_20251205.R`
- **Figure:** `FigureS9.R`

## Shared Helper Scripts

- `helpers/sim_methods.R`: shared simulation utilities used across the scenario scripts.
- `helpers/sim_nulloutcome.R`: helper for generating simulated negative-control outcomes used in the falsification tests.
- `helpers/create_nlmr_summary2025.R`: function to run uncorrected and corrected doubly-ranked method.
- `helpers/custom_summary.R`: custom summary helper.
