# Cleaning Log

Written by `01_cleaning.ipynb`. Rerun the notebook to update it.

| Dataset | Rule | Reason | Records before | Records after | Records affected |
|---|---|---|---|---|---|
| D1 | Drop duplicate ids | Each patient should appear once | 70,000 | 70,000 | 0 |
| D1 | Keep women only | Track focus on women | 70,000 | 45,530 | 24,470 |
| D1 | Require two positive BP readings | The analysis requires both systolic and diastolic BP in mmHg. | 45,530 | 45,515 | 15 |
| D1 | Exclude BP readings of 1,000 mmHg or higher | Inspection showed repeated four-digit entries and five-digit extremes. | 45,515 | 44,961 | 554 |
| D1 | Require systolic BP above diastolic BP | Equal or reversed pairs cannot be reliably corrected from this file. | 44,961 | 44,764 | 197 |
| D1 | Flag remaining BP extremes for review | Flagged 37 rows; retained all pending review. Review thresholds in mmHg: systolic <70 or >250, diastolic <40 or >150. | 44,764 | 44,764 | 0 |
| D1 | Retain rows regardless of height and weight | The primary models use age, blood pressure and cholesterol. | 44,764 | 44,764 | 0 |
| D1 | Set BP eligibility for the primary analysis | 44,727 eligible. Retained 37 flagged records in the CSV for sensitivity analysis. Primary BP bounds are 70–250 mmHg systolic and 40–150 mmHg diastolic, inclusive, with systolic above diastolic. | 44,764 | 44,764 | 0 |
| D1 | Derive age, BMI, cholesterol labels and glucose labels | Age in years = age in days / 365.25. BMI = weight in kg / height in metres squared. Cholesterol and glucose codes 1, 2, 3 map to normal, above normal, well above normal. Height and weight extremes were retained; BMI is unused in the primary analysis. | 44,764 | 44,764 | 0 |
