# SN-attributable-risk-childhood-cancer

## Overview

This repository provides R code, data processing, and statistical analysis pipelines to estimate the **Contributions of cancer treatment and genetic predisposition to risk of subsequent neoplasms in long-term survivors of childhood cancer from our publication: https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(25)00157-3/fulltext.
The focus is on quantifying the proportion of SNs attributable to treatments (chemotherapy, radiotherapy), polygenic risk scores (PRS), and lifestyle factors, including extensive stratification by age, sex, and genetic ancestry.  
**Statistical modeling and effect modification analyses** are central to the repository. 

---

## Background

Childhood cancer survivors face increased risk for subsequent neoplasms (SNs) due to a mixture of treatment exposures, inherited genetic risk, and lifestyle factors. This repository aims to rigorously estimate the **attributable fraction**—the proportion of SNs in the population that can be ascribed to specific exposures—using robust statistical methods. Analyses are performed separately for different exposures and are stratified by key modifiers (age, sex, genetic ancestry).

---

## Statistical Pipeline

### 1. Data Preparation

- **Merging Data:**  
  - Phenotype (`PHENO.ANY_SN`) and subneoplasm event data (`subneoplasms.sas7bdat`) are merged by subject ID (`sjlid`).
  - Covariates include treatment details, genetic ancestry, lifestyle factors, and demographic info.
  - Data cleaning removes duplicates, irrelevant columns (e.g., MRN), and handles missing values.

- **Variable Construction:**  
  - Exposures are categorized (e.g., dose bins for radiotherapy, chemotherapy).
  - Genetic ancestry is classified using admixture proportions (EUR: >0.8, AFR: >0.6).

- **Saving Analysis-Ready Data:**  
  - Cleaned and merged data are saved as `.RData` or `.rds` files for subsequent modeling.

### 2. Attributable Fraction Estimation

#### Core Concept

The **attributable fraction (AF)** quantifies the portion of SNs that would be prevented if a specific exposure were eliminated. Mathematically,  
AF = (Incidence in total population - Incidence in unexposed) / Incidence in total population

#### Implementation

- **Counterfactual Modeling:**  
  For each exposure (e.g., radiation, PRS, unfavorable lifestyle), the model predicts risk for each individual under "exposed" and "unexposed" scenarios.

- **Risk Calculation:**  
  - For each subgroup (e.g., EUR ancestry, AFR ancestry, age <35, ≥35), sum predicted SNs under observed and counterfactual scenarios.
  - AF for each group:  
    ```
    AF = (N_total - N_counterfactual_unexposed) / N_total
    ```
    Where:
      - N_total = predicted SNs under observed exposure
      - N_counterfactual_unexposed = predicted SNs if the exposure was absent

- **Exposure Combinations:**  
  - Combined AFs calculated for multiple exposures (e.g., joint effect of treatment + PRS + lifestyle).

- **Code Example:**  
  In scripts like `4.model_fit_Any_SN.R`, AFs are calculated and bound by group:
  ```r
  N_no_tx = sum(dat_all$pred_no_tx[dat_all$admixture == "EUR"], na.rm = TRUE)
  af_by_tx.EUR = (N_all.EUR - N_no_tx) / N_all.EUR
  af_by_tx.EUR.save = rbind(af_by_tx.EUR.save, af_by_tx.EUR)
  ```

### 3. Effect Modification Analyses

- **Purpose:**  
  Investigate whether the effect of one exposure (e.g., treatment) on SN risk is modified by another variable (e.g., ancestry, sex, age).

- **Implementation:**  
  - Stratified risk modeling: AFs are computed within subgroups (e.g., males vs. females, EUR vs. AFR ancestry).
  - Interaction terms may be included in regression models to formally test effect modification.

### 4. Subgroup Analyses

- **By Age:**  
  AFs are estimated for age groups such as <35 and ≥35 years old.

- **By Sex:**  
  Separate AFs for males and females.

- **By Genetic Ancestry:**  
  AFs for EUR and AFR ancestry groups, using admixture classification.

- **By SN Type:**  
  Specific SNs (e.g., NMSC, breast cancer, thyroid cancer, meningioma, sarcoma) are analyzed individually.

### 5. Uncertainty Quantification

- **Bootstrap Confidence Intervals:**  
  - AFs are computed across bootstrap samples to estimate variability.
  - 95% confidence intervals are calculated using percentiles of bootstrap AF estimates.

- **Code Example:**  
  ```r
  af_by_tx.EUR = paste0(round(af_by_tx.EUR.save[1,], 3), "_", 
                        round(quantile(af_by_tx.EUR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                        round(quantile(af_by_tx.EUR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  ```

### 6. Visualization

- **Plotting Results:**  
  - AF estimates and their confidence intervals are visualized using `ggplot2`.
  - Plots include:
    - AF by exposure and subgroup
    - Combined AFs for multiple exposures
    - SN-type-specific AFs

- **Saving Figures:**  
  - High-resolution plots (PDF, TIFF) are saved for publication-ready output.
  - Example:
    ```r
    ggsave(plot_name, p, width = 14.6, height = 10, device = cairo_pdf)
    ```

---

## Key Scripts

- **2.analysis_ready_samples_treatments.R:**  
  Data merging, cleaning, and preparation for statistical modeling.

- **4.model_fit_Any_SN.R:**  
  Main AF calculations; subgroup and effect modification analyses; bootstrap CIs.

- **5.model_fit_Any_SN_Effect_Modification.R:**  
  Detailed modeling of effect modification.

- **6.AF_plot.R:**  
  Advanced visualization and plot generation.

---

## Usage

1. **Prepare Data:**  
   - Ensure all phenotype, treatment, and admixture files are available.
   - Edit file paths as needed.

2. **Run Analysis:**  
   - Execute scripts in sequence:  
     `2.analysis_ready_samples_treatments.R` → `4.model_fit_Any_SN.R` → `5.model_fit_Any_SN_Effect_Modification.R` → `6.AF_plot.R`

3. **Interpret Results:**  
   - Review AF estimates and confidence intervals for exposures and subgroups.
   - Use plots for visualization and publication.

---

## Dependencies

- R (>= 3.6.0)
- `haven`, `expss`, `ggplot2`, `dplyr`
- Additional packages as referenced in scripts

---

## Citing

If you use this repository, cite the original publication or acknowledge the code and authors.

Neupane, A., et al. Contributions of cancer treatment and genetic predisposition to risk of subsequent neoplasms in long-term survivors of childhood cancer: a report from the St Jude Lifetime Cohort and the Childhood Cancer Survivor Study. The Lancet Oncology 26, 806-816 (2025).