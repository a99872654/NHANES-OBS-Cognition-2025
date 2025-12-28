# NHANES-OBS-Cognition-2025
# NHANES 2011-2014 Analysis: Oxidative Balance Score (OBS) and Cognitive Performance

This repository contains the R analytic code for the research study investigating the dose-response relationship between the **Oxidative Balance Score (OBS)** and **cognitive function** in older adults (≥60 years), using data from the **National Health and Nutrition Examination Survey (NHANES) 2011–2014**.

## 📌 Project Overview

This study utilizes complex survey-weighted models to analyze the association between a composite OBS (comprising 16 dietary and 4 lifestyle factors) and global/domain-specific cognitive performance.

**Key Analytical Features:**
* **Data Cleaning:** Rigorous exclusion criteria (age, incomplete cognition/dietary data, extreme energy intake, stroke history).
* **OBS Construction:** Calculation of OBS based on dietary and lifestyle components.
* **Statistical Modeling:** Survey-weighted linear regression (using the `survey` package) to account for NHANES complex sampling design (strata, PSU, weights).
* **Non-linearity Analysis:** Restricted Cubic Spline (RCS) analysis to test for linear vs. non-linear dose-response relationships.
* **Visualization:** Generation of flowcharts, forest plots, and trend plots.

## 🛠️ Prerequisites

To run this code, you need **R** installed. The analysis relies on the following R packages. Please ensure they are installed before running the script:

```r
install.packages(c("tidyverse", "haven", "survey", "Hmisc", "broom", "DiagrammeR", "tableone"))
