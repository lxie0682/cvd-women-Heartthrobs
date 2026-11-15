"""Reproduces every number and figure in docs/article.md.

Run from anywhere inside the repository:   python track-02/team/src/article_analysis.py
It reads the raw Dataset 1 file, applies the cleaning rules in docs/cleaning_log.md, uses the 02_analysis.ipynb split
(20% held out, random_state 42, stratified by age group and CVD label) and writes figures/, outputs/ and docs/results_summary.md.
"""

# Setup: run this cell first. It finds the track folder and the data.
import json, re, sys, subprocess, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import rankdata, chi2
from statsmodels.stats.proportion import proportion_confint
from sklearn.model_selection import train_test_split, StratifiedKFold

try:
    display
except NameError:            # running as a plain script
    display = print
warnings.filterwarnings("ignore", category=FutureWarning)

PUBLIC_REPO = "https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026.git"
IN_COLAB = "google.colab" in sys.modules

def find_track_dir():
    here = Path.cwd().resolve()
    for folder in [here, *here.parents]:
        if (folder / "track-02" / "data").is_dir():
            return folder / "track-02"
        if folder.name == "track-02" and (folder / "data").is_dir():
            return folder
    if IN_COLAB:
        dest = Path("/content/openbio-datathon-2026")
        if not dest.exists():
            subprocess.run(["git", "clone", "--depth", "1", PUBLIC_REPO, str(dest)], check=True)
        return dest / "track-02"
    return None

TRACK = find_track_dir()
# The team folder is the folder above notebooks/, wherever the submission folder has been copied.
HERE = Path.cwd().resolve()
TEAM = HERE.parent if HERE.name == "notebooks" else (TRACK / "team" if TRACK else HERE)
CLEAN, FIGURES, OUTPUTS, DOCS = TEAM / "data_clean", TEAM / "figures", TEAM / "outputs", TEAM / "docs"
for folder in (CLEAN, FIGURES, OUTPUTS, DOCS):
    folder.mkdir(parents=True, exist_ok=True)

print("Track folder:", TRACK)
print("Team folder: ", TEAM)


QUESTION = ("Among records coded as women in the OpenBio Cardiovascular Disease Dataset, how much does adding "
            "blood pressure and cholesterol improve discrimination of recorded CVD beyond age alone, and how does this "
            "improvement differ across predefined age groups?")
DATASET = "d1"
OUTCOME = "cardio"
MODEL = "logistic"

FORMULA_A = "cardio ~ age_years"                                    # Model A: age alone (baseline)
FORMULA_B = "cardio ~ age_years + ap_hi + ap_lo + C(cholesterol)"   # Model B: age + systolic BP + diastolic BP + cholesterol (primary)

AGE_BINS = [0, 45, 55, np.inf]
AGE_LABELS = ["Under 45", "45 to 54", "55 and older"]
ALL = "All women"
GROUPS = AGE_LABELS + [ALL]

SEED = 42             # same seed and age-group-by-outcome stratification as 02_analysis.ipynb
TEST_SIZE = 0.20      # 80% development, 20% held-out test
N_BOOT = 2000         # bootstrap repetitions for the primary analysis
N_BOOT_SENS = 500     # bootstrap repetitions for each sensitivity analysis

print(QUESTION)


RAW_NAME = "cardiovascular_dataset_OpenBio_Final.csv"

def load_raw():
    for folder in [TRACK / "data" if TRACK else None, Path.cwd(), Path("/content")]:
        if folder and (folder / RAW_NAME).is_file():
            return pd.read_csv(folder / RAW_NAME, sep=";")
    raise FileNotFoundError(f"{RAW_NAME} not found. Put it in track-02/data/.")

def build_cohort(raw, ap_hi=(70, 250), ap_lo=(40, 150), repair=False, drop_repeats=False, require_body=False):
    """Apply the cleaning rules from 01_cleaning.ipynb to the raw file. The defaults are the primary rules;
    the other arguments switch on the alternative rules used in the sensitivity analyses."""
    d = raw.drop_duplicates(subset="id").copy()
    d["age_years"] = d["age"] / 365.25
    d = d[d["gender"] == 1].copy()
    if repair:   # repair obvious typing errors instead of removing them
        hi, lo = d["ap_hi"].astype(float).abs(), d["ap_lo"].astype(float).abs()
        hi = hi.where(~hi.between(7, 25), hi * 10)            # 12 -> 120
        hi = hi.where(~hi.between(700, 2500), hi / 10)        # 1400 -> 140
        lo = lo.where(~lo.between(400, 2000), lo / 10)        # 1000 -> 100
        lo = lo.where(~lo.between(4, 15), lo * 10)            # 8 -> 80
        swap = (hi < lo) & hi.between(*ap_lo) & lo.between(*ap_hi)   # the two numbers typed in the wrong boxes
        hi, lo = hi.where(~swap, lo), lo.where(~swap, d["ap_hi"].astype(float).abs().where(~swap, hi))
        d["ap_hi"], d["ap_lo"] = hi, lo
    d = d[d["ap_hi"].between(*ap_hi) & d["ap_lo"].between(*ap_lo) & (d["ap_hi"] > d["ap_lo"])]
    d = d[d["cardio"].isin([0, 1]) & d["cholesterol"].isin([1, 2, 3]) & d["age_years"].between(18, 100)]
    measure_cols = [c for c in raw.columns if c != "id"]
    d["repeat_profile"] = d.duplicated(subset=measure_cols).astype(int)
    d["body_plausible"] = (d["height"].between(120, 210) & d["weight"].between(30, 200)).astype(int)
    if drop_repeats:
        d = d[d["repeat_profile"] == 0]
    if require_body:
        d = d[d["body_plausible"] == 1]
    d["age_group"] = pd.cut(d["age_years"], bins=AGE_BINS, labels=AGE_LABELS, right=False)
    return d.reset_index(drop=True)

raw = load_raw()
rebuilt = build_cohort(raw)

# The primary cohort is rebuilt from the raw file with the rules in the cleaning log (systolic 70 to 250, diastolic 40 to 150,
# systolic above diastolic). If the team's d1_clean.csv is present, confirm that its eligible records are the same women.
df = rebuilt.copy()
clean_path = CLEAN / "d1_clean.csv"
if clean_path.exists():
    saved = pd.read_csv(clean_path)
    if "bp_primary_eligible" in saved.columns:
        saved = saved[saved["bp_primary_eligible"].astype(str).str.lower().isin(["1", "true", "1.0"])]
    same = set(saved["id"]) == set(df["id"])
    print(f"d1_clean.csv found: its {len(saved):,} primary-eligible records " + ("match" if same else "DO NOT match") + f" the {len(df):,} rebuilt records.")
else:
    print(f"Primary cohort rebuilt from the raw file: {len(df):,} records.")

df["age_group"] = pd.Categorical(df["age_group"].astype(str), categories=AGE_LABELS, ordered=True)

# Integrity checks on the analysis cohort
assert (df["gender"] == 1).all() and df["id"].is_unique
assert df["ap_hi"].between(70, 250).all() and df["ap_lo"].between(40, 150).all() and (df["ap_hi"] > df["ap_lo"]).all()
assert not df[["age_years", "ap_hi", "cholesterol", "cardio"]].isna().any().any()
print(f"Analysis cohort: {len(df):,} records coded as women, {int(df['cardio'].sum()):,} labeled CVD-positive.")


# One visual identity for every figure.
INK, MUTED, SURFACE, GRID = "#14213D", "#5B6475", "#FBF8F1", "#E4DFD3"
CORAL, TEAL = "#D9623F", "#0097A7"                       # Model A (age only), Model B (age + BP + cholesterol)
AGE_COLORS = dict(zip(AGE_LABELS, ["#8FAEE0", "#4C78BF", "#1E3A75"]))   # one hue, light to dark, for the ordered age groups
AGE_MARKERS = dict(zip(AGE_LABELS, ["o", "s", "D"]))     # shape as well as color, so color is never the only cue
NO_CVD, YES_CVD = "#B9C0CC", "#1E3A75"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED, "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.axisbelow": True, "font.size": 12, "axes.titlesize": 15, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.labelsize": 12, "legend.frameon": False, "figure.dpi": 110, "savefig.dpi": 200, "font.family": "DejaVu Sans",
})

def save(fig, name):
    fig.savefig(FIGURES / name, dpi=200, bbox_inches="tight")
    print("Saved", FIGURES / name)

def wilson(k, n):
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return lo, hi

def footer(fig, text, y=-0.02):
    fig.text(0.01, y, text, fontsize=9.5, color=MUTED, ha="left", va="top")


def describe_group(g):
    q = g["age_years"].quantile([0.25, 0.5, 0.75])
    chol = g["cholesterol"].value_counts(normalize=True).reindex([1, 2, 3]).fillna(0) * 100
    return pd.Series({
        "Records": len(g), "CVD-positive": int(g["cardio"].sum()), "CVD-negative": int((1 - g["cardio"]).sum()),
        "CVD-positive %": round(100 * g["cardio"].mean(), 1),
        "Age, median (IQR)": f"{q[0.5]:.1f} ({q[0.25]:.1f} to {q[0.75]:.1f})",
        "Systolic BP, median (IQR)": f"{g['ap_hi'].median():.0f} ({g['ap_hi'].quantile(.25):.0f} to {g['ap_hi'].quantile(.75):.0f})",
        "Cholesterol normal %": round(chol[1], 1), "above normal %": round(chol[2], 1), "well above normal %": round(chol[3], 1)})

table1 = pd.concat([df.groupby("age_group", observed=True).apply(describe_group), describe_group(df).to_frame(ALL).T])
table1.to_csv(OUTPUTS / "table1_sample.csv")
display(table1)


comp = df.groupby("age_group", observed=True)["cardio"].agg(pos="sum", n="size")
comp["neg"] = comp["n"] - comp["pos"]

fig, ax = plt.subplots(figsize=(9, 5.2))
x = np.arange(len(comp))
ax.bar(x, comp["neg"], width=0.55, color=NO_CVD, edgecolor=SURFACE, linewidth=2, label="No CVD label")
ax.bar(x, comp["pos"], width=0.55, bottom=comp["neg"], color=YES_CVD, edgecolor=SURFACE, linewidth=2, label="CVD label")
for i, (g, r) in enumerate(comp.iterrows()):
    ax.text(i, r["neg"] / 2, f"{int(r['neg']):,}", ha="center", va="center", color=INK, fontsize=11)
    ax.text(i, r["neg"] + r["pos"] / 2, f"{int(r['pos']):,}", ha="center", va="center", color="white", fontsize=11)
    ax.text(i, r["n"] + 350, f"{int(r['n']):,} women", ha="center", va="bottom", fontsize=11, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(comp.index)
ax.set_xlabel("Age group"); ax.set_ylabel("Number of exam records")
ax.set_ylim(0, comp["n"].max() * 1.13)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
ax.grid(axis="x", visible=False)
ax.legend(loc="upper left")
ax.set_title("Every Age Group Has Thousands of Records With and Without a CVD Label", pad=14)
footer(fig, f"Records coded as women in Dataset 1 after cleaning (n = {len(df):,}).")
save(fig, "01_sample_composition.png"); plt.show()


band_edges = AGE_BINS
band_labels = AGE_LABELS
df["age_band5"] = pd.cut(df["age_years"], bins=band_edges, labels=band_labels, right=False)
by_band = df.groupby("age_band5", observed=True)["cardio"].agg(pos="sum", n="size")
by_band["share"] = by_band["pos"] / by_band["n"]
by_band["lo"], by_band["hi"] = wilson(by_band["pos"], by_band["n"])
display((by_band.assign(share=lambda t: (100 * t["share"]).round(1), lo=lambda t: (100 * t["lo"]).round(1), hi=lambda t: (100 * t["hi"]).round(1))))

fig, ax = plt.subplots(figsize=(8.2, 5.2))
x = np.arange(len(by_band))
ax.set_xlim(-0.5, len(by_band) - 0.5)
ax.errorbar(x, 100 * by_band["share"], yerr=[100 * (by_band["share"] - by_band["lo"]), 100 * (by_band["hi"] - by_band["share"])],
            fmt="o-", color=YES_CVD, markersize=9, linewidth=2, capsize=4, markeredgecolor=SURFACE, markeredgewidth=1.5)
for i, r in enumerate(by_band.itertuples()):
    ax.text(i, 100 * r.hi + 2.2, f"{100 * r.share:.1f}%", ha="center", fontsize=12, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"{b}\n(n = {n:,})" for b, n in zip(by_band.index, by_band["n"])], fontsize=10.5)
ax.set_ylim(0, 85); ax.set_ylabel("Records labeled with CVD (%)"); ax.set_xlabel("Age at the exam (years)")
ax.grid(axis="x", visible=False)
ax.set_title("Recorded CVD Was More Common in Older Age Groups", pad=14)
footer(fig, "Points are the share of records labeled CVD-positive in each age group. Bars are 95% Wilson confidence intervals.\n"
            "This is the share within this dataset, not the rate of CVD among women in any population.", y=-0.06)
save(fig, "02_cvd_by_age.png"); plt.show()


SBP_EDGES = [0, 120, 130, 140, 160, np.inf]
SBP_LABELS = ["Under 120", "120 to 129", "130 to 139", "140 to 159", "160 or more"]
CHOL_LABELS = {1: "Normal", 2: "Above normal", 3: "Well above normal"}
df["sbp_band"] = pd.cut(df["ap_hi"], bins=SBP_EDGES, labels=SBP_LABELS, right=False)
df["chol_label"] = pd.Categorical(df["cholesterol"].map(CHOL_LABELS), categories=list(CHOL_LABELS.values()), ordered=True)

def share_table(by):
    t = df.groupby(["age_group", by], observed=True)["cardio"].agg(pos="sum", n="size").reset_index()
    t["share"] = t["pos"] / t["n"]
    t["lo"], t["hi"] = wilson(t["pos"], t["n"])
    return t

shares = {"Cholesterol category": share_table("chol_label"), "Systolic blood pressure (mmHg)": share_table("sbp_band")}
for name, t in shares.items():
    print(name); display(t.assign(share=lambda d: (100 * d["share"]).round(1)).pivot(index=t.columns[1], columns="age_group", values="share"))
    display(t.pivot(index=t.columns[1], columns="age_group", values="n"))

fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), sharey=True)
for ax, (name, t) in zip(axes, shares.items()):
    by = t.columns[1]
    cats = list(t[by].cat.categories)
    for j, g in enumerate(AGE_LABELS):
        s = t[t["age_group"] == g].set_index(by).reindex(cats)
        xs = np.arange(len(cats)) + (j - 1) * 0.12
        ax.errorbar(xs, 100 * s["share"], yerr=[100 * (s["share"] - s["lo"]), 100 * (s["hi"] - s["share"])],
                    fmt=AGE_MARKERS[g] + "-", color=AGE_COLORS[g], markersize=8, linewidth=1.8, capsize=3,
                    markeredgecolor=SURFACE, markeredgewidth=1.2, label=g)
    ax.set_xticks(np.arange(len(cats))); ax.set_xticklabels(cats, fontsize=10.5)
    ax.set_xlabel(name); ax.set_ylim(0, 100); ax.grid(axis="x", visible=False)
axes[0].set_ylabel("Records labeled with CVD (%)")
axes[0].legend(title="Age group", loc="upper left")
fig.suptitle("At Every Age, Higher Cholesterol and Higher Blood Pressure Go With More CVD Labels",
             x=0.01, ha="left", fontsize=15, fontweight="bold", color=INK)
footer(fig, "Points are the share of records labeled CVD-positive. Bars are 95% Wilson confidence intervals. Blood pressure bands are for display only; "
            "the models use the measured value.\nDescriptive associations in exam records. They do not show that either measurement causes CVD.", y=-0.03)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save(fig, "03_measurement_associations.png"); plt.show()


ols = smf.ols("ap_hi ~ age_years", data=df).fit()
slope10, slope_ci = 10 * ols.params["age_years"], 10 * ols.conf_int().loc["age_years"]
print(f"Systolic BP rises {slope10:.1f} mmHg per 10 years of age (95% CI {slope_ci[0]:.1f} to {slope_ci[1]:.1f}); R² = {ols.rsquared:.3f}")

grid = pd.DataFrame({"age_years": np.linspace(df["age_years"].min(), df["age_years"].max(), 100)})
pred = ols.get_prediction(grid).summary_frame(alpha=0.05)
by_year = df.assign(age_int=np.floor(df["age_years"])).groupby("age_int")["ap_hi"].agg(["mean", "size"])
by_year = by_year[by_year["size"] >= 50]

fig, ax = plt.subplots(figsize=(9, 5.2))
ax.scatter(by_year.index + 0.5, by_year["mean"], s=46, color=AGE_COLORS["45 to 54"], edgecolor=SURFACE, linewidth=1.2, zorder=3,
           label="Mean for each year of age")
ax.plot(grid["age_years"], pred["mean"], color=INK, linewidth=2, label="Fitted line")
ax.fill_between(grid["age_years"], pred["mean_ci_lower"], pred["mean_ci_upper"], color=INK, alpha=0.18, label="95% confidence band")
ax.set_xlabel("Age at the exam (years)"); ax.set_ylabel("Systolic blood pressure (mmHg)")
ax.legend(loc="upper left")
ax.set_title("Older Women in the Data Have Higher Blood Pressure", pad=14)
footer(fig, f"Linear regression of systolic blood pressure on age, n = {len(df):,}. Slope: {slope10:.1f} mmHg per 10 years "
            f"(95% CI {slope_ci[0]:.1f} to {slope_ci[1]:.1f}). R² = {ols.rsquared:.2f}:\nage explains about {100 * ols.rsquared:.0f}% of the variation, "
            "so most of the difference in blood pressure between women is not about age.", y=-0.04)
save(fig, "07_bp_by_age.png"); plt.show()


dev, test = train_test_split(df, test_size=TEST_SIZE, random_state=SEED, stratify=df["age_group"].astype(str) + "_" + df["cardio"].astype(str))
dev, test = dev.reset_index(drop=True), test.reset_index(drop=True)
assert set(dev["id"]).isdisjoint(test["id"])
print(f"Development: {len(dev):,} records. Held-out test: {len(test):,} records.")
display(test.groupby("age_group", observed=True)["cardio"].agg(test_records="size", cvd_positive="sum"))

def fit(formula, data):
    return smf.logit(formula, data=data).fit(disp=False)

model_a, model_b = fit(FORMULA_A, dev), fit(FORMULA_B, dev)
test["score_a"], test["score_b"] = model_a.predict(test), model_b.predict(test)
dev["score_a"], dev["score_b"] = model_a.predict(dev), model_b.predict(dev)
print(model_b.summary().tables[1])


def auc(y, s):
    """ROC-AUC as the share of positive-negative pairs ranked correctly (ties count half)."""
    y = np.asarray(y); r = rankdata(s)
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def group_masks(data, group_col="age_group", labels=AGE_LABELS):
    masks = {g: (data[group_col].astype(str) == g).values for g in labels}
    masks[ALL] = np.ones(len(data), dtype=bool)
    return masks

def paired_bootstrap(data, score_cols, n_boot, seed, group_col="age_group", labels=AGE_LABELS):
    """Stratified paired bootstrap. Records are resampled with replacement inside each age group x label cell, so every
    replicate keeps the original number of positive and negative records in each age group. Both models' scores stay
    attached to each resampled record. Returns {score_col: {group: array of n_boot AUCs}}."""
    rng = np.random.default_rng(seed)
    y = data["cardio"].values
    g = data[group_col].astype(str).values
    cells = {(lab, k): np.flatnonzero((g == lab) & (y == k)) for lab in labels for k in (0, 1)}
    scores = {c: data[c].values for c in score_cols}
    out = {c: {grp: np.empty(n_boot) for grp in labels + [ALL]} for c in score_cols}
    for b in range(n_boot):
        idx = {key: rng.choice(ix, size=len(ix), replace=True) for key, ix in cells.items()}
        per_group = {lab: np.concatenate([idx[(lab, 0)], idx[(lab, 1)]]) for lab in labels}
        per_group[ALL] = np.concatenate(list(per_group.values()))
        for grp, ix in per_group.items():
            yy = y[ix]
            for c in score_cols:
                out[c][grp][b] = auc(yy, scores[c][ix])
    return out

def ci(a):
    return np.percentile(a, [2.5, 97.5])

def auc_table(data, boot, base="score_a", new="score_b", group_col="age_group", labels=AGE_LABELS):
    rows = []
    for grp, m in group_masks(data, group_col, labels).items():
        y = data.loc[m, "cardio"].values
        a, b = auc(y, data.loc[m, base]), auc(y, data.loc[m, new])
        d = boot[new][grp] - boot[base][grp]
        rows.append({"age_group": grp, "n_test": int(m.sum()), "n_positive": int(y.sum()), "n_negative": int((1 - y).sum()),
                     "auc_age": a, "auc_age_lo": ci(boot[base][grp])[0], "auc_age_hi": ci(boot[base][grp])[1],
                     "auc_expanded": b, "auc_expanded_lo": ci(boot[new][grp])[0], "auc_expanded_hi": ci(boot[new][grp])[1],
                     "delta_auc": b - a, "ci_lower": ci(d)[0], "ci_upper": ci(d)[1]})
    return pd.DataFrame(rows)

boot = paired_bootstrap(test, ["score_a", "score_b"], N_BOOT, SEED)
results = auc_table(test, boot)
results.to_csv(OUTPUTS / "model_results.csv", index=False)
display(results.round(3))


delta = {g: boot["score_b"][g] - boot["score_a"][g] for g in GROUPS}
pairs = [("Under 45", "45 to 54"), ("Under 45", "55 and older"), ("45 to 54", "55 and older")]
res_ix = results.set_index("age_group")
contrasts = []
for g1, g2 in pairs:
    for what, arr, point in [("Gain from adding measurements (delta AUC)", delta[g1] - delta[g2], res_ix.loc[g1, "delta_auc"] - res_ix.loc[g2, "delta_auc"]),
                             ("Expanded model AUC", boot["score_b"][g1] - boot["score_b"][g2], res_ix.loc[g1, "auc_expanded"] - res_ix.loc[g2, "auc_expanded"])]:
        lo, hi = ci(arr)
        contrasts.append({"quantity": what, "contrast": f"{g1} minus {g2}", "difference": point, "ci_lower": lo, "ci_upper": hi,
                          "interval_excludes_zero": bool(lo > 0 or hi < 0)})
contrasts = pd.DataFrame(contrasts)
contrasts.to_csv(OUTPUTS / "age_group_contrasts.csv", index=False)
display(contrasts.round(3))


def headline_figure(res, title, subtitle):
    order = [ALL] + AGE_LABELS[::-1]          # plotted bottom to top
    r = res.set_index("age_group").loc[order]
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    ypos = np.arange(len(order), dtype=float); ypos[0] -= 0.35     # set "All women" slightly apart
    for y, (g, row) in zip(ypos, r.iterrows()):
        ax.plot([row["auc_age"], row["auc_expanded"]], [y, y], color=MUTED, linewidth=2.2, zorder=1, solid_capstyle="round")
        ax.errorbar(row["auc_age"], y, xerr=[[row["auc_age"] - row["auc_age_lo"]], [row["auc_age_hi"] - row["auc_age"]]],
                    fmt="o", color=CORAL, markersize=12, capsize=4, elinewidth=1.6, markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
        ax.errorbar(row["auc_expanded"], y, xerr=[[row["auc_expanded"] - row["auc_expanded_lo"]], [row["auc_expanded_hi"] - row["auc_expanded"]]],
                    fmt="D", color=TEAL, markersize=11, capsize=4, elinewidth=1.6, markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=3)
        ax.text(0.872, y + 0.10, f"+{row['delta_auc']:.2f}", fontsize=15, fontweight="bold", va="center", ha="left")
        ax.text(0.872, y - 0.20, f"({row['ci_lower']:.2f} to {row['ci_upper']:.2f})", fontsize=10.5, color=MUTED, va="center", ha="left")
    ax.axvline(0.5, color=MUTED, linestyle=(0, (4, 4)), linewidth=1.2)
    ax.text(0.494, ypos[-1] + 0.62, "coin flip", ha="right", va="bottom", fontsize=10, color=MUTED)
    ax.text(0.872, ypos[-1] + 0.62, "Gain (95% CI)", ha="left", va="bottom", fontsize=10.5, color=MUTED)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{g}\n{int(row['n_positive']):,} with / {int(row['n_negative']):,} without CVD label" for g, row in r.iterrows()], fontsize=11)
    for lab in ax.get_yticklabels():
        lab.set_linespacing(1.5)
    ax.set_xlim(0.46, 0.95); ax.set_ylim(ypos[0] - 0.6, ypos[-1] + 0.95)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8])
    ax.set_xlabel("How often the model ranks a record with CVD above a record without (AUC)")
    ax.grid(axis="y", visible=False)
    handles = [plt.Line2D([], [], marker="o", color=CORAL, linestyle="", markersize=11, markeredgecolor=SURFACE, label="Age alone"),
               plt.Line2D([], [], marker="D", color=TEAL, linestyle="", markersize=10, markeredgecolor=SURFACE, label="Age + blood pressure (both numbers) + cholesterol")]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, -0.30), ncol=2, fontsize=11)
    fig.suptitle(title, x=0.01, y=1.05, ha="left", fontsize=17, fontweight="bold", color=INK)
    fig.text(0.01, 0.975, subtitle, ha="left", va="top", fontsize=11.5, color=MUTED, linespacing=1.5)
    footer(fig, f"Held-out test set of {int(res.set_index('age_group').loc[ALL, 'n_test']):,} records coded as women (OpenBio Dataset 1). Horizontal bars and brackets are 95% confidence intervals "
                f"from {N_BOOT:,} paired bootstrap resamples.\nThe models rank records by a CVD label recorded at the exam. They do not predict future heart disease. Not medical advice or a diagnostic tool.", y=-0.13)
    return fig

r_ = results.set_index("age_group")
top = r_.loc[AGE_LABELS, "delta_auc"].idxmax()
HEADLINE_TITLE = "Blood Pressure and Cholesterol Tell Us Less About Recorded CVD in Women 55 and Older"
HEADLINE_SUB = (f"Adding the two measurements to age improves how well a model separates records with and without a CVD label in every age group:\n"
                f"by {r_.loc['Under 45', 'delta_auc']:.2f} under 45, {r_.loc['45 to 54', 'delta_auc']:.2f} at 45 to 54 and {r_.loc['55 and older', 'delta_auc']:.2f} at 55 and older.")
fig = headline_figure(results, HEADLINE_TITLE, HEADLINE_SUB)
save(fig, "headline.png"); save(fig, "04_incremental_auc.png"); plt.show()


def logit_fn(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))

def calibration_stats(y, p):
    lp = logit_fn(np.asarray(p))
    slope = sm.Logit(y, sm.add_constant(lp)).fit(disp=False).params[1]
    intercept = sm.GLM(y, np.ones(len(y)), family=sm.families.Binomial(), offset=lp).fit().params[0]
    return {"brier": float(np.mean((p - y) ** 2)), "calibration_slope": float(slope), "calibration_intercept": float(intercept),
            "mean_score": float(np.mean(p)), "observed_share": float(np.mean(y))}

cal_rows = []
for grp, m in group_masks(test).items():
    for name, col in [("Age alone", "score_a"), ("Age + BP + cholesterol", "score_b")]:
        cal_rows.append({"age_group": grp, "model": name, **calibration_stats(test.loc[m, "cardio"].values, test.loc[m, col].values)})
calibration = pd.DataFrame(cal_rows)
calibration.to_csv(OUTPUTS / "calibration.csv", index=False)
display(calibration.round(3))

def calibration_curve(y, p, bins=10):
    q = pd.qcut(p, bins, duplicates="drop")
    t = pd.DataFrame({"y": y, "p": p, "q": q}).groupby("q", observed=True).agg(mean_p=("p", "mean"), pos=("y", "sum"), n=("y", "size"))
    t["obs"] = t["pos"] / t["n"]; t["lo"], t["hi"] = wilson(t["pos"], t["n"])
    return t

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8), sharey=True)
ax = axes[0]
for name, col, color, mk in [("Age alone", "score_a", CORAL, "o"), ("Age + blood pressure + cholesterol", "score_b", TEAL, "D")]:
    t = calibration_curve(test["cardio"].values, test[col].values)
    ax.errorbar(100 * t["mean_p"], 100 * t["obs"], yerr=[100 * (t["obs"] - t["lo"]), 100 * (t["hi"] - t["obs"])], fmt=mk + "-", color=color,
                markersize=8, linewidth=1.8, capsize=3, markeredgecolor=SURFACE, markeredgewidth=1.2, label=name)
ax.set_title("Both Models, All Women", fontsize=13)
ax = axes[1]
for g in AGE_LABELS:
    m = (test["age_group"].astype(str) == g).values
    t = calibration_curve(test.loc[m, "cardio"].values, test.loc[m, "score_b"].values, bins=6)
    ax.errorbar(100 * t["mean_p"], 100 * t["obs"], yerr=[100 * (t["obs"] - t["lo"]), 100 * (t["hi"] - t["obs"])], fmt=AGE_MARKERS[g] + "-",
                color=AGE_COLORS[g], markersize=8, linewidth=1.8, capsize=3, markeredgecolor=SURFACE, markeredgewidth=1.2, label=g)
ax.set_title("Expanded Model, by Age Group", fontsize=13)
for ax in axes:
    ax.plot([0, 100], [0, 100], color=MUTED, linestyle=(0, (4, 4)), linewidth=1.2)
    ax.text(84, 90, "perfect\nagreement", fontsize=9.5, color=MUTED, ha="center", va="bottom", rotation=0)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_xlabel("Model score (%)"); ax.legend(loc="upper left", fontsize=10.5)
axes[0].set_ylabel("Records actually labeled with CVD (%)")
fig.suptitle("Scores Match Reality Overall, but Are Too Cautious for Younger Women and Too Bold for Older Women", x=0.01, ha="left", fontsize=15, fontweight="bold", color=INK)
footer(fig, "Held-out test records are sorted by model score and split into equal-sized bins. Points on the dashed line mean that, for example, about 30% of the records "
            "scored near 30%\nwere labeled with CVD. A line steeper than the dashed one means scores are bunched too close together. Bars are 95% Wilson confidence intervals. Agreement with a recorded label is not validation for predicting future disease.", y=-0.02)
fig.tight_layout(rect=(0, 0, 1, 0.95))
save(fig, "05_calibration.png"); plt.show()


SENS = []   # one row per analysis x age group

def run_analysis(label, kind, data, formula_a=FORMULA_A, formula_b=FORMULA_B, n_boot=N_BOOT_SENS, note=""):
    """Split with the same seed and stratification, fit both models on the development part, and bootstrap the test part."""
    d, t = train_test_split(data, test_size=TEST_SIZE, random_state=SEED, stratify=data["age_group"].astype(str) + "_" + data["cardio"].astype(str))
    d, t = d.reset_index(drop=True), t.reset_index(drop=True)
    t["score_a"], t["score_b"] = fit(formula_a, d).predict(t), fit(formula_b, d).predict(t)
    tab = auc_table(t, paired_bootstrap(t, ["score_a", "score_b"], n_boot, SEED))
    tab.insert(0, "analysis", label); tab.insert(1, "kind", kind); tab["n_cohort"] = len(data); tab["note"] = note
    SENS.append(tab)
    return tab

primary_row = results.copy(); primary_row.insert(0, "analysis", "Primary analysis"); primary_row.insert(1, "kind", "primary")
primary_row["n_cohort"] = len(df); primary_row["note"] = "Rules, models and age groups as prespecified"
SENS.append(primary_row)

# 1. Flexible (spline) terms for age and blood pressure, in both models
run_analysis("Flexible curves for age and blood pressure", "model form", df,
             "cardio ~ cr(age_years, df=4)", "cardio ~ cr(age_years, df=4) + cr(ap_hi, df=4) + cr(ap_lo, df=4) + C(cholesterol)",
             note="Natural cubic splines, 4 df each; the age-only model gets the same flexible age term")

# 2. Restore the 37 flagged blood pressure records. Existing women keep their original development/test assignment.
full = build_cohort(raw, ap_hi=(1, 999), ap_lo=(1, 999))
flagged = full[~full["id"].isin(df["id"])].copy()
N_FLAGGED = len(flagged)
f_dev, f_test = train_test_split(flagged, test_size=TEST_SIZE, random_state=SEED)
d2 = pd.concat([dev[full.columns.intersection(dev.columns)], f_dev[full.columns.intersection(dev.columns)]], ignore_index=True)
t2 = pd.concat([test[full.columns.intersection(test.columns)], f_test[full.columns.intersection(test.columns)]], ignore_index=True)
t2["age_group"] = pd.Categorical(t2["age_group"].astype(str), categories=AGE_LABELS, ordered=True)
t2["score_a"], t2["score_b"] = fit(FORMULA_A, d2).predict(t2), fit(FORMULA_B, d2).predict(t2)
tab = auc_table(t2, paired_bootstrap(t2, ["score_a", "score_b"], N_BOOT_SENS, SEED))
tab.insert(0, "analysis", f"{N_FLAGGED} flagged blood pressure records restored"); tab.insert(1, "kind", "cleaning rule")
tab["n_cohort"] = len(full); tab["note"] = "Existing women keep their original assignment"
SENS.append(tab)

variants = {
    "Stricter blood pressure range": dict(ap_hi=(90, 200), ap_lo=(50, 120)),
    "Typing errors repaired, not removed": dict(repair=True),
}
for label, kwargs in variants.items():
    cohort = build_cohort(raw, **kwargs)
    note = ", ".join(f"{k}={v}" for k, v in kwargs.items())
    run_analysis(label, "cleaning rule", cohort, note=f"{note}; {len(cohort):,} records ({len(cohort) - len(df):+,} vs primary)")

# 3. Systolic pressure only (no diastolic) in the expanded model
run_analysis("Expanded model without diastolic pressure", "model form", df, FORMULA_A, "cardio ~ age_years + ap_hi + C(cholesterol)")

sens = pd.concat(SENS, ignore_index=True)
display(sens[sens["age_group"] != ALL].pivot(index="analysis", columns="age_group", values="delta_auc").loc[
    [s["analysis"].iloc[0] for s in SENS], AGE_LABELS].round(3))
display(sens.groupby("analysis", sort=False)["n_cohort"].first().to_frame("records in cohort"))


# 3. Repeated cross-validation: does the result depend on which records landed in the test set?
N_REPEATS, N_FOLDS = 20, 5
strata_all = (df["age_group"].astype(str) + "_" + df["cardio"].astype(str)).values
cv_rows = []
for rep in range(N_REPEATS):
    oof = pd.DataFrame({"cardio": df["cardio"].values, "age_group": df["age_group"].astype(str).values, "a": np.nan, "b": np.nan})
    for tr, te in StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED + rep).split(df, strata_all):
        oof.loc[te, "a"] = fit(FORMULA_A, df.iloc[tr]).predict(df.iloc[te]).values
        oof.loc[te, "b"] = fit(FORMULA_B, df.iloc[tr]).predict(df.iloc[te]).values
    for grp, m in group_masks(oof).items():
        y = oof.loc[m, "cardio"].values
        cv_rows.append({"repeat": rep, "age_group": grp, "auc_age": auc(y, oof.loc[m, "a"]), "auc_expanded": auc(y, oof.loc[m, "b"])})
cv = pd.DataFrame(cv_rows); cv["delta_auc"] = cv["auc_expanded"] - cv["auc_age"]
cv_summary = cv.groupby("age_group", sort=False).agg(delta_mean=("delta_auc", "mean"), delta_min=("delta_auc", "min"), delta_max=("delta_auc", "max"),
                                                      auc_age_mean=("auc_age", "mean"), auc_expanded_mean=("auc_expanded", "mean")).loc[GROUPS]
display(cv_summary.round(3))
cv_row = pd.DataFrame({"analysis": f"Cross-validation, all records ({N_REPEATS} x {N_FOLDS}-fold)", "kind": "resampling", "age_group": GROUPS,
                       "auc_age": cv_summary["auc_age_mean"].values, "auc_expanded": cv_summary["auc_expanded_mean"].values,
                       "delta_auc": cv_summary["delta_mean"].values, "ci_lower": cv_summary["delta_min"].values, "ci_upper": cv_summary["delta_max"].values,
                       "n_cohort": len(df), "note": "Mean over repeats; the range shown is min to max across repeats, not a confidence interval"})
sens = pd.concat([sens, cv_row], ignore_index=True)


# 4. Do the age-group cut points matter? Five-year bands and equal-sized thirds.
test["band5"] = pd.cut(test["age_years"], bins=[0, 45, 50, 55, 60, np.inf], labels=["Under 45", "45 to 49", "50 to 54", "55 to 59", "60 to 64"], right=False).astype(str)
labels5 = ["Under 45", "45 to 49", "50 to 54", "55 to 59", "60 to 64"]
boot5 = paired_bootstrap(test, ["score_a", "score_b"], N_BOOT_SENS, SEED, group_col="band5", labels=labels5)
bands5 = auc_table(test, boot5, group_col="band5", labels=labels5)

cuts = dev["age_years"].quantile([1 / 3, 2 / 3]).values          # thirds defined on the development set
labels3 = [f"Youngest third (under {cuts[0]:.1f})", f"Middle third ({cuts[0]:.1f} to {cuts[1]:.1f})", f"Oldest third ({cuts[1]:.1f} and over)"]
test["third"] = pd.cut(test["age_years"], bins=[0, cuts[0], cuts[1], np.inf], labels=labels3, right=False).astype(str)
boot3 = paired_bootstrap(test, ["score_a", "score_b"], N_BOOT_SENS, SEED, group_col="third", labels=labels3)
thirds = auc_table(test, boot3, group_col="third", labels=labels3)

alt_bands = pd.concat([bands5.assign(scheme="Five-year bands"), thirds.assign(scheme="Equal-sized thirds")])
alt_bands = alt_bands[alt_bands["age_group"] != ALL]
alt_bands.to_csv(OUTPUTS / "alternative_age_bands.csv", index=False)
display(alt_bands[["scheme", "age_group", "n_test", "n_positive", "auc_age", "auc_expanded", "delta_auc", "ci_lower", "ci_upper"]].round(3))


sens.to_csv(OUTPUTS / "sensitivity_analyses.csv", index=False)
order = list(dict.fromkeys(sens["analysis"]))[::-1]
prim = results.set_index("age_group")

fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.6), sharey=True, sharex=True)
for ax, g in zip(axes, AGE_LABELS):
    s = sens[sens["age_group"] == g].set_index("analysis").loc[order]
    ax.axvspan(prim.loc[g, "ci_lower"], prim.loc[g, "ci_upper"], color=TEAL, alpha=0.12, linewidth=0)
    ax.axvline(prim.loc[g, "delta_auc"], color=TEAL, linewidth=1.2, linestyle=(0, (4, 3)))
    for i, (name, row) in enumerate(s.iterrows()):
        is_p = name == "Primary analysis"
        ax.errorbar(row["delta_auc"], i, xerr=[[row["delta_auc"] - row["ci_lower"]], [row["ci_upper"] - row["delta_auc"]]],
                    fmt="D" if is_p else "o", color=TEAL if is_p else INK, markersize=10 if is_p else 7.5, capsize=3.5, elinewidth=1.5,
                    markeredgecolor=SURFACE, markeredgewidth=1.2)
    ax.set_title(g, fontsize=13); ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=10.5)
    ax.grid(axis="y", visible=False); ax.set_xlabel("Gain in AUC from adding\nblood pressure and cholesterol")
fig.suptitle("The Result Holds Under Every Alternative We Tried", x=0.01, ha="left", fontsize=15, fontweight="bold", color=INK)
footer(fig, "Each row repeats the whole analysis with one decision changed. Points are the gain in held-out AUC; bars are 95% bootstrap confidence intervals "
            "(for cross-validation, the range across 20 repeats).\nThe shaded band is the primary analysis's confidence interval.", y=-0.03)
fig.tight_layout(rect=(0, 0, 1, 0.95))
save(fig, "08_sensitivity.png"); plt.show()


EXPLORATORY = {
    "Age + blood pressure": "cardio ~ age_years + ap_hi + ap_lo",
    "Age + cholesterol": "cardio ~ age_years + C(cholesterol)",
    "Age + blood pressure + cholesterol (primary)": FORMULA_B,
    "Age + everything recorded": "cardio ~ age_years + ap_hi + ap_lo + C(cholesterol) + C(gluc) + smoke + alco + active + weight",
}
expl_cols = {}
for i, (name, f) in enumerate(EXPLORATORY.items()):
    col = f"score_x{i}"; test[col] = fit(f, dev).predict(test); expl_cols[name] = col
boot_x = paired_bootstrap(test, ["score_a"] + list(expl_cols.values()), N_BOOT_SENS, SEED)
expl = []
for name, col in expl_cols.items():
    t = auc_table(test, boot_x, new=col); t.insert(0, "model", name); expl.append(t)
expl = pd.concat(expl, ignore_index=True)
expl.insert(1, "analysis_type", np.where(expl["model"].str.contains("primary"), "primary", "exploratory"))
expl.to_csv(OUTPUTS / "exploratory_models.csv", index=False)
display(expl.pivot(index="model", columns="age_group", values="auc_expanded").loc[list(EXPLORATORY), GROUPS].round(3))
display(expl.pivot(index="model", columns="age_group", values="delta_auc").loc[list(EXPLORATORY), GROUPS].round(3))


cut = test["score_b"].quantile(0.25)
pos = test[test["cardio"] == 1].copy()
pos["low_scored"] = pos["score_b"] <= cut
low = pos.groupby("age_group", observed=True).agg(cvd_positive=("id", "size"), low_scored=("low_scored", "sum"))
low["low_scored_%"] = (100 * low["low_scored"] / low["cvd_positive"]).round(1)
lo_ci = [wilson(k, n) for k, n in zip(low["low_scored"], low["cvd_positive"])]
low["ci_lower_%"] = [round(100 * c[0], 1) for c in lo_ci]; low["ci_upper_%"] = [round(100 * c[1], 1) for c in lo_ci]
print(f"Lowest-quartile cut-off for the expanded model's score: {cut:.3f}")
display(low)

profile = pos.groupby("low_scored").agg(records=("id", "size"), median_age=("age_years", "median"), median_systolic=("ap_hi", "median"),
                                        normal_cholesterol_pct=("cholesterol", lambda s: 100 * (s == 1).mean()),
                                        systolic_under_130_pct=("ap_hi", lambda s: 100 * (s < 130).mean()))
profile.index = profile.index.map({True: "CVD-positive, low-scored", False: "CVD-positive, other"})
display(profile.round(1))

quiet = (df["ap_hi"] < 130) & (df["cholesterol"] == 1)
df["quiet"] = quiet
unseen = df[df["cardio"] == 1].groupby("age_group", observed=True)["quiet"].agg(k="sum", n="size")
unseen["share_%"] = (100 * unseen["k"] / unseen["n"]).round(1)
uc = [wilson(k, n) for k, n in zip(unseen["k"], unseen["n"])]
unseen["ci_lower_%"] = [round(100 * c[0], 1) for c in uc]; unseen["ci_upper_%"] = [round(100 * c[1], 1) for c in uc]
print("CVD-positive records with systolic BP under 130 mmHg and normal cholesterol, whole cleaned cohort:")
display(unseen)
low.to_csv(OUTPUTS / "low_scored_positives.csv"); unseen.to_csv(OUTPUTS / "unremarkable_positives.csv")


df["age_per5"], df["sbp_per10"], df["dbp_per10"] = df["age_years"] / 5, df["ap_hi"] / 10, df["ap_lo"] / 10
INFER = "cardio ~ age_per5 + sbp_per10 + dbp_per10 + C(cholesterol)"

def or_table(res, label):
    t = pd.DataFrame({"OR": np.exp(res.params), "ci_lower": np.exp(res.conf_int()[0]), "ci_upper": np.exp(res.conf_int()[1]), "p_value": res.pvalues})
    t = t.drop(index="Intercept").rename(index={"age_per5": "Age, per 5 years", "sbp_per10": "Systolic BP, per 10 mmHg", "dbp_per10": "Diastolic BP, per 10 mmHg",
                                                 "C(cholesterol)[T.2]": "Cholesterol above normal (vs normal)",
                                                 "C(cholesterol)[T.3]": "Cholesterol well above normal (vs normal)"})
    t.insert(0, "cohort", label); t.insert(1, "n", int(res.nobs))
    return t.rename_axis("term").reset_index()

overall_fit = fit(INFER, df)
ors = [or_table(overall_fit, ALL)] + [or_table(fit(INFER, df[df["age_group"] == g]), g) for g in AGE_LABELS]
ors = pd.concat(ors, ignore_index=True)
ors.to_csv(OUTPUTS / "odds_ratios.csv", index=False)
display(ors.round(3))

# Does the association differ by age group? Likelihood-ratio test of age group x measurement interactions.
base_fit = fit("cardio ~ age_per5 + C(age_group) + sbp_per10 + dbp_per10 + C(cholesterol)", df)
int_sbp = fit("cardio ~ age_per5 + C(age_group) * sbp_per10 + dbp_per10 + C(cholesterol)", df)
int_chol = fit("cardio ~ age_per5 + C(age_group) * C(cholesterol) + sbp_per10 + dbp_per10", df)
def lr_test(small, big):
    stat = 2 * (big.llf - small.llf); dof = int(big.df_model - small.df_model)
    return stat, dof, chi2.sf(stat, dof)
lr_sbp, lr_chol = lr_test(base_fit, int_sbp), lr_test(base_fit, int_chol)
print(f"Age group x systolic BP:  LR chi-square = {lr_sbp[0]:.1f} on {lr_sbp[1]} df, p = {lr_sbp[2]:.2g}")
print(f"Age group x cholesterol:  LR chi-square = {lr_chol[0]:.1f} on {lr_chol[1]} df, p = {lr_chol[2]:.2g}")


terms = ["Age, per 5 years", "Systolic BP, per 10 mmHg", "Diastolic BP, per 10 mmHg", "Cholesterol above normal (vs normal)", "Cholesterol well above normal (vs normal)"]
fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6), gridspec_kw={"width_ratios": [1, 1.25]})

ax = axes[0]
o = ors[ors["cohort"] == ALL].set_index("term").loc[terms[::-1]]
ax.errorbar(o["OR"], range(len(o)), xerr=[o["OR"] - o["ci_lower"], o["ci_upper"] - o["OR"]], fmt="D", color=TEAL, markersize=10, capsize=4,
            elinewidth=1.6, markeredgecolor=SURFACE, markeredgewidth=1.2)
for i, r in enumerate(o.itertuples()):
    ax.text(r.ci_upper * 1.04, i, f"{r.OR:.2f} ({r.ci_lower:.2f} to {r.ci_upper:.2f})", va="center", fontsize=10.5)
ax.set_yticks(range(len(o))); ax.set_yticklabels(["Age,\nper 5 years", "Systolic BP,\nper 10 mmHg", "Diastolic BP,\nper 10 mmHg", "Cholesterol\nabove normal", "Cholesterol well\nabove normal"][::-1])
ax.set_title(f"All Women (n = {len(df):,})", fontsize=13)

ax = axes[1]
sub_terms = terms[1:][::-1]
for j, g in enumerate(AGE_LABELS):
    o = ors[ors["cohort"] == g].set_index("term").loc[sub_terms]
    ys = np.arange(len(sub_terms)) + (1 - j) * 0.22
    ax.errorbar(o["OR"], ys, xerr=[o["OR"] - o["ci_lower"], o["ci_upper"] - o["OR"]], fmt=AGE_MARKERS[g], color=AGE_COLORS[g], markersize=8.5,
                capsize=3.5, elinewidth=1.5, markeredgecolor=SURFACE, markeredgewidth=1.2, label=g)
ax.set_yticks(range(len(sub_terms))); ax.set_yticklabels(["Systolic BP,\nper 10 mmHg", "Diastolic BP,\nper 10 mmHg", "Cholesterol\nabove normal", "Cholesterol well\nabove normal"][::-1])
ax.set_title("By Age Group (Exploratory)", fontsize=13); ax.legend(title="Age group", loc="upper right", fontsize=10.5)
for ax in axes:
    ax.set_xscale("log"); ax.axvline(1, color=MUTED, linestyle=(0, (4, 4)), linewidth=1.2)
    ax.set_xticks([1, 1.5, 2, 3, 4, 6]); ax.set_xticklabels(["1", "1.5", "2", "3", "4", "6"]); ax.minorticks_off()
    ax.set_xlabel("Odds ratio for a recorded CVD label (log scale)"); ax.grid(axis="y", visible=False)
axes[0].set_xlim(0.9, 5.2); axes[1].set_xlim(0.9, 11)
fig.suptitle("Blood Pressure and Cholesterol Are Tied to a CVD Label Most Strongly in Younger Women",
             x=0.01, ha="left", fontsize=15, fontweight="bold", color=INK)
footer(fig, "Logistic regression on the cleaned cohort; each odds ratio is adjusted for the other variables shown. Bars are 95% confidence intervals. An odds ratio of 1 (dashed line) means no association.\n"
            "Right panel: the same model fitted separately within each age group, also adjusted for age. Associations in exam records, not causal effects.", y=-0.03)
fig.tight_layout(rect=(0, 0, 1, 0.95))
save(fig, "06_regression_forest.png"); plt.show()


r = results.set_index("age_group"); c = contrasts.set_index(["quantity", "contrast"]); cal = calibration.set_index(["age_group", "model"])
o = ors.set_index(["cohort", "term"])
def f3(x): return f"{x:.3f}"
def rng(lo, hi, d=3): return f"{lo:.{d}f} to {hi:.{d}f}"

lines = ["# Results Summary", "",
         "Written by `notebooks/02_analysis.ipynb`. Rerun the notebook to update it. Every number comes from the code.", "",
         "> This project is a student data analysis. It is not medical advice, a diagnostic tool or clinical guidance.", "",
         "## Question", "", QUESTION, "",
         "## Sample", "",
         f"- Cleaned cohort: **{len(df):,}** records coded as women, {int(df['cardio'].sum()):,} labeled CVD-positive ({100 * df['cardio'].mean():.1f}%).",
         f"- Development set: {len(dev):,} records (80%). Held-out test set: {len(test):,} records (20%). Random split stratified by age group and CVD label, seed {SEED} (the split in 02_analysis.ipynb).", "",
         "## Models", "",
         f"- Model A (baseline): `{FORMULA_A}`", f"- Model B (expanded, primary): `{FORMULA_B}`",
         "- Both are unpenalized logistic regressions fitted once on the development set and evaluated on the same held-out records.", "",
         "## Primary result: gain in held-out AUC from adding blood pressure and cholesterol", "",
         f"95% confidence intervals from {N_BOOT:,} paired bootstrap resamples, stratified by age group and CVD label.", "",
         "| Age group | Test records | CVD-positive | CVD-negative | AUC, age alone (95% CI) | AUC, expanded (95% CI) | Gain (95% CI) |", "|---|---:|---:|---:|---|---|---|"]
for g in GROUPS:
    x = r.loc[g]
    lines.append(f"| {g} | {int(x['n_test']):,} | {int(x['n_positive']):,} | {int(x['n_negative']):,} | {f3(x['auc_age'])} ({rng(x['auc_age_lo'], x['auc_age_hi'])}) | "
                 f"{f3(x['auc_expanded'])} ({rng(x['auc_expanded_lo'], x['auc_expanded_hi'])}) | **+{f3(x['delta_auc'])}** ({rng(x['ci_lower'], x['ci_upper'])}) |")
lines += ["", "## Direct comparisons between age groups", "",
          "A difference between two groups is supported only when the interval for the difference itself excludes zero.", "",
          "| Quantity | Contrast | Difference (95% CI) | Interval excludes zero |", "|---|---|---|---|"]
for (q, k), x in c.iterrows():
    lines.append(f"| {q} | {k} | {x['difference']:+.3f} ({rng(x['ci_lower'], x['ci_upper'])}) | {'yes' if x['interval_excludes_zero'] else 'no'} |")
lines += ["", "## Calibration on the held-out test set", "",
          "A slope of 1 and an intercept of 0 mean scores agree with observed frequencies. The Brier score is the mean squared gap between score and label (lower is better; 0.25 is what a constant 50% score earns).", "",
          "| Age group | Model | Brier score | Calibration slope | Calibration intercept |", "|---|---|---:|---:|---:|"]
for (g, m), x in cal.iterrows():
    lines.append(f"| {g} | {m} | {x['brier']:.3f} | {x['calibration_slope']:.2f} | {x['calibration_intercept']:+.2f} |")
lines += ["", "## Sensitivity analyses: gain in AUC by age group", "",
          "Each row repeats the whole pipeline (split, fit, bootstrap) with one decision changed. Cleaning variants change which records exist, so they also re-draw the split.", "",
          "| Analysis | Records | " + " | ".join(AGE_LABELS) + " |", "|---|---:|---|---|---|"]
for name in dict.fromkeys(sens["analysis"]):
    s_ = sens[sens["analysis"] == name].set_index("age_group")
    cells = [f"+{s_.loc[g, 'delta_auc']:.3f} ({rng(s_.loc[g, 'ci_lower'], s_.loc[g, 'ci_upper'])})" for g in AGE_LABELS]
    lines.append(f"| {name} | {int(s_['n_cohort'].iloc[0]):,} | " + " | ".join(cells) + " |")
lines += ["", "For the cross-validation row the bracket is the range across repeats, not a confidence interval.", "",
          "## Alternative age bands (same models, same test set)", "", "| Scheme | Band | Test records | AUC, age alone | AUC, expanded | Gain (95% CI) |", "|---|---|---:|---:|---:|---|"]
for x in alt_bands.itertuples():
    lines.append(f"| {x.scheme} | {x.age_group} | {x.n_test:,} | {x.auc_age:.3f} | {x.auc_expanded:.3f} | +{x.delta_auc:.3f} ({rng(x.ci_lower, x.ci_upper)}) |")
lines += ["", "## Exploratory: which measurement carries the gain?", "", "| Model | Type | " + " | ".join(GROUPS) + " |", "|---|---|" + "---|" * len(GROUPS)]
for name in EXPLORATORY:
    e_ = expl[expl["model"] == name].set_index("age_group")
    lines.append(f"| {name} | {e_['analysis_type'].iloc[0]} | " + " | ".join(f"{e_.loc[g, 'auc_expanded']:.3f} (+{e_.loc[g, 'delta_auc']:.3f})" for g in GROUPS) + " |")
lines += ["", "Cells show held-out AUC and, in brackets, the gain over age alone.", "",
          "## Adjusted associations (odds ratios, 95% CI)", "", "| Cohort | Term | Odds ratio (95% CI) |", "|---|---|---|"]
for (coh, term), x in o.iterrows():
    lines.append(f"| {coh} | {term} | {x['OR']:.2f} ({rng(x['ci_lower'], x['ci_upper'], 2)}) |")
lines += ["", f"Likelihood-ratio tests for a difference across age groups: systolic blood pressure chi-square = {lr_sbp[0]:.1f} on {lr_sbp[1]} df (p = {lr_sbp[2]:.1g}); "
          f"cholesterol chi-square = {lr_chol[0]:.1f} on {lr_chol[1]} df (p = {lr_chol[2]:.1g}). The by-age-group odds ratios are exploratory.", "",
          "## Exploratory: CVD-positive records the expanded model scores low", "",
          f"Low-scored means a score in the lowest quarter of the held-out cohort (at or below {cut:.3f}). This is a description of model errors, not a diagnostic threshold, and these records are not 'missed diagnoses'.", "",
          "| Age group | CVD-positive test records | Low-scored | Share (95% CI) |", "|---|---:|---:|---|"]
for g, x in low.iterrows():
    lines.append(f"| {g} | {int(x['cvd_positive']):,} | {int(x['low_scored']):,} | {x['low_scored_%']:.1f}% ({x['ci_lower_%']:.1f} to {x['ci_upper_%']:.1f}) |")
lines += ["", "Share of all CVD-positive records in the cleaned cohort with systolic blood pressure under 130 mmHg and normal cholesterol:", "",
          "| Age group | CVD-positive records | Under 130 mmHg and normal cholesterol | Share (95% CI) |", "|---|---:|---:|---|"]
for g, x in unseen.iterrows():
    lines.append(f"| {g} | {int(x['n']):,} | {int(x['k']):,} | {x['share_%']:.1f}% ({x['ci_lower_%']:.1f} to {x['ci_upper_%']:.1f}) |")
lines += ["", "## Context", "", f"Systolic blood pressure rises {slope10:.1f} mmHg per 10 years of age (95% CI {slope_ci[0]:.1f} to {slope_ci[1]:.1f}); R² = {ols.rsquared:.3f}.", "",
          "## Limits on interpretation", "",
          "- The outcome is a CVD label recorded at the same exam as the measurements. Nothing here predicts future disease.",
          "- The source does not document how the label was assigned. If it reflects diagnoses that include high blood pressure itself, blood pressure is partly measuring the outcome.",
          "- Associations are not causal effects. Age groups are age bands, not menopause stages.",
          "- The bootstrap intervals describe uncertainty from the test sample. They do not cover how the records were collected, which the source does not describe.",
          "- Sex is recorded as one binary code. The cohort is a convenience sample, not a sample of any population of women.", ""]
(DOCS / "results_summary.md").write_text("\n".join(lines), encoding="utf-8")
print("Wrote", DOCS / "results_summary.md")

key = {"n_clean": len(df), "n_test": len(test), "results": results.round(4).to_dict(orient="records"), "contrasts": contrasts.round(4).to_dict(orient="records")}
(OUTPUTS / "key_numbers.json").write_text(json.dumps(key, indent=1), encoding="utf-8")


u, m, ol, al = (r.loc[g] for g in GROUPS)
cu = c.loc[("Gain from adding measurements (delta AUC)", "Under 45 minus 55 and older")]
cm = c.loc[("Gain from adding measurements (delta AUC)", "Under 45 minus 45 to 54")]
print(f"""ANSWER TO OUR QUESTION

Among {len(df):,} records coded as women, adding blood pressure and cholesterol to age improved how well a
logistic model separated records with and without a CVD label in every age group, on {len(test):,} records the model never saw.

  Under 45:     AUC {u['auc_age']:.2f} -> {u['auc_expanded']:.2f}, gain +{u['delta_auc']:.2f} (95% CI {u['ci_lower']:.2f} to {u['ci_upper']:.2f})
  45 to 54:     AUC {m['auc_age']:.2f} -> {m['auc_expanded']:.2f}, gain +{m['delta_auc']:.2f} (95% CI {m['ci_lower']:.2f} to {m['ci_upper']:.2f})
  55 and older: AUC {ol['auc_age']:.2f} -> {ol['auc_expanded']:.2f}, gain +{ol['delta_auc']:.2f} (95% CI {ol['ci_lower']:.2f} to {ol['ci_upper']:.2f})

The gain was larger under 45 than at 55 and older (difference {cu['difference']:+.2f}, 95% CI {cu['ci_lower']:.2f} to {cu['ci_upper']:.2f}).
The difference between under 45 and 45 to 54 ({cm['difference']:+.2f}, 95% CI {cm['ci_lower']:.2f} to {cm['ci_upper']:.2f}) includes zero, so we do not claim one.

These are associations with a label recorded at the same exam, in one dataset. They are not predictions of future
heart disease, not causal effects, and not medical advice.""")


