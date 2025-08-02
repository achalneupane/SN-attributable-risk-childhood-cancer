# SN-Attributable-Risk-Childhood-Cancer

## Overview

This repository contains R scripts and workflows for analyzing the **attributable risk of subsequent neoplasms (SNs)** in long-term survivors of childhood cancer, as described in our [Lancet Oncology publication (2025)](https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045%2825%2900157-3/fulltext).

We quantify the contribution of cancer treatments (chemotherapy, radiotherapy), polygenic risk scores (PRS), and lifestyle factors to SN risk—stratified by age, sex, and genetic ancestry—using robust statistical modeling and effect modification analyses.

---

## Background

Childhood cancer survivors are at elevated risk for SNs due to a combination of treatment exposures, inherited genetic predisposition, and modifiable lifestyle factors. This repository aims to estimate the **attributable fraction (AF)**—the proportion of SNs that could be avoided if a specific exposure were removed—through counterfactual modeling and survival analysis techniques.

---

## Repository Structure & Workflow

### 1. Data Preparation

**Scripts:**

* `1.data_preparation_with_lifestyle_factors.R`
* `2.analysis_ready_samples_treatments.R`

**Steps:**

* Merge phenotype, treatment, and lifestyle data by subject ID (`sjlid`).
* Clean data: remove duplicates, handle missing values, and drop irrelevant columns.
* Construct variables (e.g., categorize chemotherapy/radiotherapy doses, infer genetic ancestry from admixture proportions).
* Save clean, analysis-ready datasets as `.RData` or `.rds`.

---

### 2. Time-to-Event Formatting for SN Analysis

**Script:**

* `3.SJLIFE_SN_Longitudinal_Data_Preparation_TimeSplit_V1.R`

**Description:**
Prepares longitudinal SN data for survival analysis by splitting follow-up into yearly intervals. It captures recurrent events, defines risk periods starting 5 years post-diagnosis, and creates spline terms for age to model time-varying hazards.

---

### 3. Attributable Fraction Estimation

**Script:**

* `4.model_fit_Any_SN.R`

This script employs a statistical approach based on Poisson regression to estimate the risk of developing subsequent neoplasms among childhood cancer survivors. The model incorporates various covariates—such as age, sex, treatment exposures (chemotherapy and radiotherapy), genetic ancestry, and lifestyle factors—to predict the expected number of neoplasms for each individual. To quantify the impact of specific exposures, the code generates counterfactual scenarios by setting exposures (e.g., radiotherapy, chemotherapy) to "None" and recalculates the predicted risks. The attributable fraction (AF) for each exposure is then calculated as the proportion of predicted cases that would be prevented if the exposure were eliminated, both in the overall population and within subgroups defined by age, sex, and ancestry. Confidence intervals for AF estimates are derived from bootstrap replicates, providing robust measures of uncertainty for each risk estimate. This detailed, subgroup-specific AF analysis enables a comprehensive understanding of the factors driving SN risk in childhood cancer survivors.

**Concept:**
AF = (Incidence\_observed – Incidence\_counterfactual\_unexposed) / Incidence\_observed

**Approach:**

* Predict individual-level SN risk under observed and counterfactual (unexposed) conditions.
* Compute AF for each subgroup (e.g., by ancestry, age, sex).
* Estimate combined AFs for multiple exposures (e.g., treatment + PRS + lifestyle).
* Store results for visualization and reporting.

---

### 4. Effect Modification Analysis

**Script:**

* `5.model_fit_Any_SN_Effect_Modification.R`

**Purpose:**
Explore whether exposure effects on SN risk differ by age, sex, or genetic ancestry.

**Method:**

* Stratified modeling and interaction terms to assess effect modification.
* Report subgroup-specific AFs.

---

### 5. Subgroup Analyses

AFs are estimated for:

* **Age:** <35 vs. ≥35 years
* **Sex:** Males vs. females
* **Genetic Ancestry:** EUR (>80% European) vs. AFR (>60% African)
* **SN Types:** Includes NMSC, breast cancer, thyroid cancer, meningioma, sarcoma, etc.

---

### 6. Uncertainty Quantification

**Bootstrap Confidence Intervals:**

* AFs are resampled across bootstrap replicates.
* 95% confidence intervals are derived from the empirical distribution.
* Stored in matrix form for downstream visualization.

---

### 7. Visualization

**Script:**

* `6.AF_plot.R`

**Features:**

* Generate high-resolution plots using `ggplot2`:

  * AF by exposure and subgroup
  * Combined and SN-type-specific AFs
* Save publication-ready figures in PDF/TIFF format.

---

## How to Use

### 1. Setup

Ensure the following inputs are available:

* Phenotype and SN event data
* Treatment and lifestyle data
* Genetic ancestry estimates

Edit file paths as needed within each script.

### 2. Run the Pipeline

Recommended script order:

```r
2.analysis_ready_samples_treatments.R  
3.SJLIFE_SN_Longitudinal_Data_Preparation_TimeSplit_V1.R  
4.model_fit_Any_SN.R  
5.model_fit_Any_SN_Effect_Modification.R  
6.AF_plot.R
```

### 3. Review Results

* Examine subgroup- and exposure-specific AF estimates with bootstrapped CIs.
* Use generated plots for visual summaries and figure panels.

---

## Dependencies

Ensure the following R packages are installed:

* `haven`
* `expss`
* `dplyr`
* `ggplot2`
* Other packages as called in scripts

Minimum R version: 3.6.0

---

## Citing

If you use this repository, cite the original publication or acknowledge the code and authors.

Neupane, A., et al. Contributions of cancer treatment and genetic predisposition to risk of subsequent neoplasms in long-term survivors of childhood cancer: a report from the St Jude Lifetime Cohort and the Childhood Cancer Survivor Study. The Lancet Oncology 26, 806-816 (2025).