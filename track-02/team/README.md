# Heartthrobs: The Unseen Heart

OpenBio 2026 Datathon, Track 2: Cardiovascular Health in Women

> This project is a student data analysis. It is not medical advice, a diagnostic tool or clinical guidance.

## Team

| Name | Role |
|---|---|
| Padmaja Krishna Singh | Team member |
| Kimmone Bartley | Team member |
| Lucy Xie | Team member |
| Vasu Jetly | Team member |

## Research Question

Among records coded as women in the OpenBio Cardiovascular Disease Dataset, how much does adding blood pressure and cholesterol improve discrimination of recorded CVD beyond age alone, and how does this improvement differ across predefined age groups (under 45, 45 to 54, 55 and older)?

## Datasets Used

- [x] Dataset 1: Cardiovascular Disease Dataset
- [ ] Dataset 2: Heart Attack Risk Prediction Dataset (synthetic)
- [ ] Dataset 3: Global CVD in Women, Tables S1 and S2

## How to Run

1. The raw files are already in `track-02/data/` (see the [data guide](../data/README.md)). On Colab the first notebook cell downloads them.
2. Install packages: `pip install -r requirements.txt`
3. Run `notebooks/01_cleaning.ipynb`, then `notebooks/02_analysis.ipynb`, top to bottom.
4. To regenerate every number, table and figure used in the article in one go (the article analysis), run `python track-02/team/src/article_analysis.py`. `notebooks/02_analysis.ipynb` points to it in its first cell. It uses the same split as `02_analysis.ipynb` (20% held out, random_state 42, stratified by age group and CVD label).

## Files

| Path | Contents |
|---|---|
| `notebooks/01_cleaning.ipynb` | Loads and cleans the data, writes the cleaning log |
| `notebooks/02_analysis.ipynb` | Summary statistics, models, checks and the headline visual |
| `src/article_analysis.py` | **Article analysis.** One script (the code behind `docs/article.md`) that reproduces every number, table and figure in the article, including the sensitivity checks, calibration and odds ratios that go beyond `02_analysis.ipynb`. Start here to check any figure in the article |
| `src/headline_bargraph.py` | Makes Figure 4 (held-out AUC bars by age group) |
| `figures/` | The six article figures (`headline.png`, `02_cvd_by_age.png`, `03_measurement_associations.png`, `headline_bargraph.png`, `06_regression_forest.png`, `08_sensitivity.png`) plus supporting graphs |
| `docs/cleaning_log.md` | Each cleaning rule and the records affected |
| `docs/article.md` | Science communication article, following the [HUDP article guide](../../docs/article-guide.md) |
| `requirements.txt` | Python packages |

## Key Findings

In 44,727 records coded as women, adding blood pressure and cholesterol to age improved held-out AUC by 0.25 (95% CI 0.21 to 0.29) under 45, 0.24 (0.21 to 0.26) at 45 to 54 and 0.15 (0.13 to 0.17) at 55 and older. The gain was larger under 45 than at 55 and older (difference 0.10, 95% CI 0.05 to 0.15), and 45 to 54 also exceeded 55 and older (0.08, 0.05 to 0.11). The gain under 45 versus 45 to 54 did not differ clearly (0.01, -0.03 to 0.06). The result held across sensitivity analyses (flexible curves, alternative cleaning rules, cross-validation, alternative age bands).

## Limitations

The outcome is a CVD label recorded at the same exam as the measurements, so nothing here predicts future disease, and the source does not document how the label was assigned. Associations are not causal effects. Sex is recorded as one binary code, age groups are not menopause stages, and the records are a convenience sample whose collection is not described. We used only Dataset 1, so the findings are not tested on other populations.

## Presentation

Interactive website (Evidence Explorer, scan or click): https://claude.ai/artifact/DuGE3TNgwHm3H6qhxzgeSB

Video link: https://youtu.be/29cp7lBAf1M
