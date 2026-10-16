# Predicting delayed graft function after DCD kidney transplantation

Code for a machine-learning model that predicts delayed graft function (DGF)
after kidney transplantation from a donation-after-circulatory-death (DCD)
donor, using UNOS STAR registry data plus donor hemodynamic measurements taken
from withdrawal of life-sustaining treatment to time of death.

Outcome: DGF, defined as dialysis in the first week after transplantation
(`FIRST_WK_DIAL`).

Cohort: DCD kidney transplants from January 1, 2015 to December 31, 2023.
Exclusions: donors missing hemodynamic data, donors aged < 18 years,
recipients aged < 18 years, multiorgan transplants, and en bloc kidney
transplants.

## Results

Cohort: 17,001 DCD kidney donors and 30,621 DCD kidney transplants.

| Model | AUROC (95% CI) |
|---|---|
| XGBoost, 5-fold CV | 0.74 (0.73-0.75) |
| Random forest, 5-fold CV | 0.73 (0.73-0.74) |
| XGBoost, test set (raw) | 0.76 (0.75-0.77) |
| XGBoost, test set (isotonic calibration) | 0.76 (0.74-0.77) |

5-fold CV values are the mean AUROC and 95% CI across held-out folds; test-set
CIs use the DeLong method. The recurrent neural network, long short-term
memory, and transformer models reported in the manuscript are not part of
this repository.

- AUPRC (XGBoost): 0.67
- Brier score: 0.20 raw vs 0.20 calibrated
- Hosmer-Lemeshow goodness of fit: p < 0.001 raw, p = 0.13 calibrated

## Data

The analysis uses UNOS Standard Transplant Analysis and Research (STAR) data.
The data cannot be redistributed, so it is not included here; it has to be
requested separately from UNOS/OPTN. Three files go in `data/`:

| File | Contents |
|---|---|
| `KIDPAN_DATA_CLEAN_MERGED.dta` | recipient / donor / transplant records |
| `DECEASED_DONOR_DCD_MEASURES.DTA` | donor hemodynamic measurements, ~1/min |
| `Diagnosis Categorization.xlsx` | `DIAG_KI` code -> category map |

## Requirements

R, plus:

```r
install.packages(c("haven", "tidyverse", "recipes", "ranger", "pROC",
                   "readxl", "ResourceSelection", "xgboost", "ggpubr", "pdp"))
```

Reading the main `.dta` file needs a few GB of RAM.

## How to run

Put the three data files in `data/`, then run the scripts from the repository
root. Every script defaults to `PROJECT_DIR <- "."`, so no path editing is
needed; change that one line only if the data lives somewhere else.

Run in this order:

```
# 1. preprocessing, cross-validation, final XGBoost, isotonic calibration
Rscript build_model.R

# 2. the expensive per-figure computations (rerun only when the model changes)
Rscript compute_feature_importance.R
Rscript compute_pdp.R
Rscript compute_permutation_importance.R

# 3. drawing (reads the CSVs from step 2, never touches the model)
Rscript plot_feature_importance.R
Rscript plot_pdp.R
Rscript plot_permutation_importance.R
```

`translate_varname.R` is not run directly; it is sourced by
`compute_feature_importance.R` and `plot_permutation_importance.R` so both
label variables identically.

## Methods notes

- Donors are split 70%/30% into training and held-out sets; the held-out set
  is then split equally into calibration (15%) and test (15%) sets. All
  transplants from one donor stay in the same partition.
- Missing values are imputed (mean for numeric variables, mode for
  categorical variables), fit on the training set and applied to the
  calibration and test sets.
- For XGBoost, categorical predictors with fewer than two levels in the
  training set are removed, and missing or unseen categories go to an
  explicit "missing" level before one-hot encoding (no intercept). The final
  model uses a binary logistic objective with class weighting and 500
  boosting rounds.
- Donor hemodynamic series are forward-filled within each donor. The summary
  features are the maximum, range, standard deviation, last recorded value,
  and linear slope over time of systolic blood pressure, diastolic blood
  pressure, and oxygen saturation.