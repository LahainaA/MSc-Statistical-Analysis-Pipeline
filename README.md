# Appendix: Statistical Analysis Pipeline (R Code)

### Master of Science in Cognitive Science
**University of Malta (UOM)**
**Student:** Laura-Viviane Engels
**Thesis Title:** Representational Momentum for Curvilinear Motion: Evidence from Sinusoidal Paths


---

## Overview
This repository serves as the digital appendix for my MSc thesis. It contains the complete **R Markdown (`.Rmd`)** source code used for the statistical analysis, data preprocessing, and visualization reported in the study.

The code is provided here to demonstrate the **analytical logic, model specifications, and reproducibility** of the results presented in Chapter 4 of the manuscript.

---

## Data Availability Statement
**Privacy & Ethics Note:**
In accordance with the University of Malta’s research ethics and data protection guidelines (GDPR), the raw participant data associated with this project is **not included** in this public repository. 

* **Data Status:** Excluded to protect participant confidentiality.
* **Code Execution:** The R Markdown files are set to `eval=FALSE` (non-executable mode). They are intended for **static inspection** of the code structure and statistical approach, rather than for active execution.

---

## Repository Structure
The analysis pipeline is organized into the following R Markdown documents:

* **`00_Demographics.Rmd`**: 
* Aggregated participant demographics (Age, Gender, Handedness). 
* Summary tables and APA-style reporting text for all three experiments. 

* **`01_Experiment_1_Analysis.Rmd`**: 
* **Preprocessing:** Data cleaning, outlier removal (>2.5 SD), and trial exclusion. 
* **Main Analysis:** Linear Mixed Effects Models (LMEM) and RM-ANOVA for the motion task. 
* **Visualization:** Generation of the primary plots and interaction effects. 

* **`02_Landmark_Exploratory.Rmd`**: 
* **Derivative Method:** Sine-wave fitting to trajectory data to identify landmark regions. 
* **Analysis:** Testing the effect of "turning points" vs. slopes on forward error. 
* **Robustness:** Sensitivity checks using different slope thresholds (10%, 15%, 20%). 

* **`03_Cross_Experiment_Comparison.Rmd`**: 
* **The "Grand Finale":** Combined analysis of Horizontal (E1), Upward (E2), and Downward (E3) experiments. * **Split-Half Validation:** Continuous moderation analysis using independent covariate baselines. 
* **Subgroup Analysis:** Comparing "High Shift" vs. "Low Shift" strategies across all motion axes.


* **`Experiment Source Code`**: 
The raw source code:** for the behavioral tasks (written in HTML/JS/CSS) is located in the `Experiment_Task` directory
---

## 🛠️ Technical Details
For examiners wishing to understand the computational environment, this analysis was built using:

* **Language:** R (Version 2025.09.2+418)
* **Key Libraries:** `tidyverse`, `lme4`, `lmerTest`, `ggplot2`, `afex`

---
