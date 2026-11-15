# Results Summary

Written by `notebooks/02_analysis.ipynb`. Rerun the notebook to update it. Every number comes from the code.

> This project is a student data analysis. It is not medical advice, a diagnostic tool or clinical guidance.

## Question

Among records coded as women in the OpenBio Cardiovascular Disease Dataset, how much does adding blood pressure and cholesterol improve discrimination of recorded CVD beyond age alone, and how does this improvement differ across predefined age groups?

## Sample

- Cleaned cohort: **44,727** records coded as women, 22,011 labeled CVD-positive (49.2%).
- Development set: 35,781 records (80%). Held-out test set: 8,946 records (20%). Random split stratified by age group and CVD label, seed 42 (the split in 02_analysis.ipynb).

## Models

- Model A (baseline): `cardio ~ age_years`
- Model B (expanded, primary): `cardio ~ age_years + ap_hi + ap_lo + C(cholesterol)`
- Both are unpenalized logistic regressions fitted once on the development set and evaluated on the same held-out records.

## Primary result: gain in held-out AUC from adding blood pressure and cholesterol

95% confidence intervals from 2,000 paired bootstrap resamples, stratified by age group and CVD label.

| Age group | Test records | CVD-positive | CVD-negative | AUC, age alone (95% CI) | AUC, expanded (95% CI) | Gain (95% CI) |
|---|---:|---:|---:|---|---|---|
| Under 45 | 1,253 | 349 | 904 | 0.585 (0.550 to 0.618) | 0.835 (0.807 to 0.862) | **+0.251** (0.211 to 0.293) |
| 45 to 54 | 3,655 | 1,591 | 2,064 | 0.549 (0.530 to 0.568) | 0.786 (0.771 to 0.800) | **+0.236** (0.215 to 0.258) |
| 55 and older | 4,038 | 2,463 | 1,575 | 0.575 (0.558 to 0.593) | 0.728 (0.713 to 0.744) | **+0.153** (0.133 to 0.173) |
| All women | 8,946 | 4,403 | 4,543 | 0.648 (0.644 to 0.653) | 0.794 (0.786 to 0.803) | **+0.146** (0.138 to 0.155) |

## Direct comparisons between age groups

A difference between two groups is supported only when the interval for the difference itself excludes zero.

| Quantity | Contrast | Difference (95% CI) | Interval excludes zero |
|---|---|---|---|
| Gain from adding measurements (delta AUC) | Under 45 minus 45 to 54 | +0.014 (-0.032 to 0.061) | no |
| Expanded model AUC | Under 45 minus 45 to 54 | +0.050 (0.019 to 0.078) | yes |
| Gain from adding measurements (delta AUC) | Under 45 minus 55 and older | +0.098 (0.053 to 0.146) | yes |
| Expanded model AUC | Under 45 minus 55 and older | +0.107 (0.075 to 0.138) | yes |
| Gain from adding measurements (delta AUC) | 45 to 54 minus 55 and older | +0.084 (0.054 to 0.114) | yes |
| Expanded model AUC | 45 to 54 minus 55 and older | +0.058 (0.037 to 0.078) | yes |

## Calibration on the held-out test set

A slope of 1 and an intercept of 0 mean scores agree with observed frequencies. The Brier score is the mean squared gap between score and label (lower is better; 0.25 is what a constant 50% score earns).

| Age group | Model | Brier score | Calibration slope | Calibration intercept |
|---|---|---:|---:|---:|
| Under 45 | Age alone | 0.198 | 2.33 | +0.01 |
| Under 45 | Age + BP + cholesterol | 0.136 | 1.67 | +0.03 |
| 45 to 54 | Age alone | 0.244 | 0.80 | -0.01 |
| 45 to 54 | Age + BP + cholesterol | 0.184 | 1.18 | -0.02 |
| 55 and older | Age alone | 0.234 | 1.23 | +0.00 |
| 55 and older | Age + BP + cholesterol | 0.204 | 0.78 | +0.01 |
| All women | Age alone | 0.233 | 1.02 | -0.00 |
| All women | Age + BP + cholesterol | 0.186 | 1.02 | +0.00 |

## Sensitivity analyses: gain in AUC by age group

Each row repeats the whole pipeline (split, fit, bootstrap) with one decision changed. Cleaning variants change which records exist, so they also re-draw the split.

| Analysis | Records | Under 45 | 45 to 54 | 55 and older |
|---|---:|---|---|---|
| Primary analysis | 44,727 | +0.251 (0.211 to 0.293) | +0.236 (0.215 to 0.258) | +0.153 (0.133 to 0.173) |
| Flexible curves for age and blood pressure | 44,727 | +0.254 (0.217 to 0.295) | +0.239 (0.218 to 0.261) | +0.154 (0.135 to 0.172) |
| 37 flagged blood pressure records restored | 44,764 | +0.250 (0.213 to 0.288) | +0.235 (0.215 to 0.256) | +0.153 (0.134 to 0.174) |
| Stricter blood pressure range | 44,619 | +0.232 (0.194 to 0.270) | +0.249 (0.228 to 0.270) | +0.150 (0.129 to 0.170) |
| Typing errors repaired, not removed | 45,477 | +0.233 (0.194 to 0.277) | +0.233 (0.211 to 0.254) | +0.166 (0.146 to 0.185) |
| Expanded model without diastolic pressure | 44,727 | +0.250 (0.213 to 0.293) | +0.237 (0.216 to 0.258) | +0.152 (0.133 to 0.172) |
| Cross-validation, all records (20 x 5-fold) | 44,727 | +0.245 (0.244 to 0.245) | +0.247 (0.246 to 0.247) | +0.153 (0.153 to 0.153) |

For the cross-validation row the bracket is the range across repeats, not a confidence interval.

## Alternative age bands (same models, same test set)

| Scheme | Band | Test records | AUC, age alone | AUC, expanded | Gain (95% CI) |
|---|---|---:|---:|---:|---|
| Five-year bands | Under 45 | 1,253 | 0.585 | 0.835 | +0.251 (0.214 to 0.293) |
| Five-year bands | 45 to 49 | 1,375 | 0.503 | 0.805 | +0.302 (0.264 to 0.336) |
| Five-year bands | 50 to 54 | 2,280 | 0.546 | 0.767 | +0.221 (0.190 to 0.247) |
| Five-year bands | 55 to 59 | 2,398 | 0.525 | 0.738 | +0.214 (0.186 to 0.241) |
| Five-year bands | 60 to 64 | 1,640 | 0.548 | 0.681 | +0.134 (0.098 to 0.169) |
| Equal-sized thirds | Youngest third (under 50.5) | 2,969 | 0.573 | 0.819 | +0.246 (0.227 to 0.270) |
| Equal-sized thirds | Middle third (50.5 to 57.3) | 2,951 | 0.548 | 0.756 | +0.208 (0.181 to 0.233) |
| Equal-sized thirds | Oldest third (57.3 and over) | 3,026 | 0.579 | 0.719 | +0.141 (0.115 to 0.165) |

## Exploratory: which measurement carries the gain?

| Model | Type | Under 45 | 45 to 54 | 55 and older | All women |
|---|---|---|---|---|---|
| Age + blood pressure | exploratory | 0.823 (+0.239) | 0.773 (+0.223) | 0.718 (+0.143) | 0.785 (+0.136) |
| Age + cholesterol | exploratory | 0.654 (+0.069) | 0.617 (+0.067) | 0.623 (+0.048) | 0.684 (+0.036) |
| Age + blood pressure + cholesterol (primary) | primary | 0.835 (+0.251) | 0.786 (+0.236) | 0.728 (+0.153) | 0.794 (+0.146) |
| Age + everything recorded | exploratory | 0.842 (+0.257) | 0.789 (+0.240) | 0.728 (+0.152) | 0.796 (+0.148) |

Cells show held-out AUC and, in brackets, the gain over age alone.

## Adjusted associations (odds ratios, 95% CI)

| Cohort | Term | Odds ratio (95% CI) |
|---|---|---|
| All women | Cholesterol above normal (vs normal) | 1.46 (1.37 to 1.55) |
| All women | Cholesterol well above normal (vs normal) | 2.69 (2.50 to 2.90) |
| All women | Age, per 5 years | 1.32 (1.30 to 1.34) |
| All women | Systolic BP, per 10 mmHg | 1.77 (1.73 to 1.81) |
| All women | Diastolic BP, per 10 mmHg | 1.13 (1.09 to 1.17) |
| Under 45 | Cholesterol above normal (vs normal) | 2.63 (2.17 to 3.19) |
| Under 45 | Cholesterol well above normal (vs normal) | 4.95 (3.67 to 6.67) |
| Under 45 | Age, per 5 years | 2.08 (1.71 to 2.55) |
| Under 45 | Systolic BP, per 10 mmHg | 2.36 (2.19 to 2.56) |
| Under 45 | Diastolic BP, per 10 mmHg | 1.13 (1.01 to 1.26) |
| 45 to 54 | Cholesterol above normal (vs normal) | 1.52 (1.38 to 1.68) |
| 45 to 54 | Cholesterol well above normal (vs normal) | 2.99 (2.64 to 3.39) |
| 45 to 54 | Age, per 5 years | 1.07 (1.01 to 1.14) |
| 45 to 54 | Systolic BP, per 10 mmHg | 1.91 (1.84 to 1.98) |
| 45 to 54 | Diastolic BP, per 10 mmHg | 1.17 (1.11 to 1.24) |
| 55 and older | Cholesterol above normal (vs normal) | 1.24 (1.13 to 1.36) |
| 55 and older | Cholesterol well above normal (vs normal) | 2.36 (2.15 to 2.58) |
| 55 and older | Age, per 5 years | 1.53 (1.45 to 1.62) |
| 55 and older | Systolic BP, per 10 mmHg | 1.57 (1.52 to 1.62) |
| 55 and older | Diastolic BP, per 10 mmHg | 1.06 (1.01 to 1.12) |

Likelihood-ratio tests for a difference across age groups: systolic blood pressure chi-square = 310.9 on 2 df (p = 3e-68); cholesterol chi-square = 84.9 on 4 df (p = 2e-17). The by-age-group odds ratios are exploratory.

## Exploratory: CVD-positive records the expanded model scores low

Low-scored means a score in the lowest quarter of the held-out cohort (at or below 0.306). This is a description of model errors, not a diagnostic threshold, and these records are not 'missed diagnoses'.

| Age group | CVD-positive test records | Low-scored | Share (95% CI) |
|---|---:|---:|---|
| Under 45 | 349 | 126 | 36.1% (31.2 to 41.3) |
| 45 to 54 | 1,591 | 191 | 12.0% (10.5 to 13.7) |
| 55 and older | 2,463 | 95 | 3.9% (3.2 to 4.7) |

Share of all CVD-positive records in the cleaned cohort with systolic blood pressure under 130 mmHg and normal cholesterol:

| Age group | CVD-positive records | Under 130 mmHg and normal cholesterol | Share (95% CI) |
|---|---:|---:|---|
| Under 45 | 1,743 | 545 | 31.3% (29.1 to 33.5) |
| 45 to 54 | 7,952 | 2,378 | 29.9% (28.9 to 30.9) |
| 55 and older | 12,316 | 3,396 | 27.6% (26.8 to 28.4) |

## Context

Systolic blood pressure rises 5.8 mmHg per 10 years of age (95% CI 5.6 to 6.1); R² = 0.054.

## Limits on interpretation

- The outcome is a CVD label recorded at the same exam as the measurements. Nothing here predicts future disease.
- The source does not document how the label was assigned. If it reflects diagnoses that include high blood pressure itself, blood pressure is partly measuring the outcome.
- Associations are not causal effects. Age groups are age bands, not menopause stages.
- The bootstrap intervals describe uncertainty from the test sample. They do not cover how the records were collected, which the source does not describe.
- Sex is recorded as one binary code. The cohort is a convenience sample, not a sample of any population of women.
