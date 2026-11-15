"""Headline bar graph: AUC of age-only vs age+BP+cholesterol model, by age group.
Reads outputs/model_results.csv (written by article_analysis.py)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

r = pd.read_csv("outputs/model_results.csv").set_index("age_group").loc[["Under 45","45 to 54","55 and older"]]
INK, MUTED, CORAL = "#11132a", "#b9b5d1", "#e4573d"
fig, ax = plt.subplots(figsize=(9, 5.6), dpi=200)
x = np.arange(3); w = 0.36
for off, col, lab, c, lo, hi in [(-w/2, "auc_age", "Age only", MUTED, "auc_age_lo", "auc_age_hi"),
                                 ( w/2, "auc_expanded", "Age + blood pressure + cholesterol", CORAL, "auc_expanded_lo", "auc_expanded_hi")]:
    v = r[col].values
    ax.bar(x+off, v, w*0.94, color=c, label=lab, zorder=3)
    ax.errorbar(x+off, v, yerr=[v-r[lo].values, r[hi].values-v], fmt="none", ecolor=INK, capsize=4, lw=1.3, zorder=4)
    for xi, vi, h in zip(x+off, v, r[hi].values):
        ax.text(xi, h+0.012, f"{vi:.3f}", ha="center", va="bottom", fontsize=11, color=INK, fontweight="bold" if c==CORAL else "normal")
for xi, (_, row) in zip(x, r.iterrows()):
    ax.text(xi, 0.10, f"gain +{row.delta_auc:.3f}\n95% CI {row.ci_lower:.3f} to {row.ci_upper:.3f}", ha="center", va="center", linespacing=1.4, fontsize=9.5, color="white", fontweight="bold", zorder=5,
            bbox=dict(boxstyle="round,pad=0.35", fc=INK, ec="none"))
ax.axhline(0.5, color=INK, ls="--", lw=1, zorder=2)
ax.set_xlim(-0.6, 2.95); ax.text(2.93, 0.51, "coin flip (0.5)", ha="right", va="bottom", fontsize=9.5, color=INK)
ax.set_xticks(x); ax.set_xticklabels([f"{g}\n(n = {int(n):,} held-out women)" for g, n in zip(r.index, r.n_test)], fontsize=10.5)
ax.set_ylim(0, 1.0); ax.set_ylabel("AUC on held-out records (higher = better)", fontsize=10.5)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e6e4ee", lw=0.8, zorder=0)
ax.legend(loc="upper right", frameon=False, fontsize=10, ncol=1, bbox_to_anchor=(1, 1.0))
fig.suptitle("Blood pressure and cholesterol add less for women 55 and older", x=0.07, ha="left", fontsize=15, fontweight="bold", color=INK, y=0.985)
fig.text(0.07, 0.915, "How well each model separates records with and without a recorded CVD label (44,727 women)", fontsize=10.5, color="#4a4a63")
fig.text(0.07, 0.012, "Error bars and gain intervals: 95% CIs from 2,000 paired bootstrap resamples.\nRecords labeled with CVD at the exam; not a prediction of future disease. Not medical advice.", fontsize=8.5, color="#4a4a63")
fig.subplots_adjust(top=0.84, bottom=0.22, left=0.09, right=0.97)
fig.savefig("figures/headline_bargraph.png"); print("saved")
