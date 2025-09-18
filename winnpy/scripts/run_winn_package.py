from pathlib import Path
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt

# === 1) Path settings ===
ROOT = Path(__file__).resolve().parent        # scripts/
DATA_DIR = ROOT.parent / "data"               # Pywinn/data

XLSX_PATH = DATA_DIR / "raetsch1_DATA_FILTERED.xlsx"
META_PATH = DATA_DIR / "per_sample_info.csv"

OUTDIR = ROOT / "results"                     # scripts/results
OUTDIR.mkdir(exist_ok=True)

# === 2) Read data (keep consistent with original script) ===
metabolite_intensities = pd.read_excel(XLSX_PATH, sheet_name="intensities1")
# Fix FutureWarning: do not use dtypes[0], use iloc or to_numpy
first_dtype = metabolite_intensities.dtypes.iloc[0] if hasattr(metabolite_intensities.dtypes, "iloc") \
              else metabolite_intensities.dtypes.to_numpy()[0]
if not np.issubdtype(first_dtype, np.number):
    metabolite_intensities = metabolite_intensities.set_index(metabolite_intensities.columns[0])
metabolite_intensities = metabolite_intensities.apply(pd.to_numeric, errors="coerce")

sample_metadata = pd.read_csv(META_PATH)

# Fix column names (consistent with your original script)
metabolite_intensities.columns = (
    metabolite_intensities.columns.astype(str)
      .str.replace("HG20.535", "HG2020.535", regex=False)
      .str.replace("HG20.709", "HG2020.709", regex=False)
)

# Reorder by sampleID
if "sampleID" not in sample_metadata.columns:
    raise ValueError("Metadata is missing 'sampleID'")
cols = sample_metadata["sampleID"].astype(str).tolist()
metabolite_intensities.columns = metabolite_intensities.columns.astype(str)
missing = set(cols) - set(metabolite_intensities.columns)
if missing:
    raise ValueError(f"These samples are not found: {sorted(missing)}")
metabolite_intensities = metabolite_intensities.loc[:, cols]

# === 3) Call WiNN ===
try:
    from winnpy import winn
except Exception:
    from winnpy.winn import winn

for need in ["plate", "order", "sample"]:
    if need not in sample_metadata.columns:
        raise ValueError(f"Metadata is missing '{need}'")

batch_codes = pd.Categorical(sample_metadata["plate"]).codes + 1
run_order = sample_metadata["order"].to_numpy()
control_cols = sample_metadata.loc[
    sample_metadata["sample"].astype(str).str.lower() == "control", "sampleID"
].astype(str).tolist() or None

params = dict(
    batch=batch_codes,
    run_order=run_order,
    control_samples=control_cols,
    parameters="fixed",
    fdr_threshold=0.05,
    median_adjustment="normalize",
    detrend_non_autocorrelated="mean",
    remove_batch_effects="combat",   # Try combat first
    test="Ljung-Box",
    lag=20,
    scale_by_batch=False,
    spline_method="conservative",
)

try:
    X = winn(data=metabolite_intensities, **params)
except Exception as e:
    print(f"[WARN] ComBat failed ({e}). Falling back to ANOVA.")
    params["remove_batch_effects"] = "anova"
    X = winn(data=metabolite_intensities, **params)

# === 4) Calculate CV and plot ===
eps = np.finfo(float).eps
if control_cols:
    X_ctrl = X.loc[:, control_cols]
    cv_control = X_ctrl.std(axis=1, ddof=1) / (X_ctrl.mean(axis=1) + eps)
else:
    cv_control = pd.Series(dtype=float, name="cv")

cv_all = X.std(axis=1, ddof=1) / (X.mean(axis=1) + eps)

cv_control_data = pd.DataFrame({"method": "winn", "group": "control", "cv": cv_control.values})
cv_all_data     = pd.DataFrame({"method": "winn", "group": "all",     "cv": cv_all.values})
cv_combined_data = pd.concat([cv_control_data, cv_all_data], ignore_index=True)

summary = (cv_combined_data
           .groupby("group")["cv"]
           .agg(mean="mean", std="std", count="count")
           .reset_index())
summary["se"] = summary["std"] / np.sqrt(summary["count"].clip(lower=1))

fig, ax = plt.subplots()
xpos = {"control": 0.9, "all": 1.1}
for _, row in summary.iterrows():
    ax.errorbar(xpos[row["group"]], row["mean"], yerr=row["se"], fmt="o", capsize=3, label=row["group"])
ax.set_xticks([1.0]); ax.set_xticklabels(["winn"], rotation=90)
ax.set_ylabel("Coefficient of Variation (CV)")
ax.set_title("CV (RSD) across sample replicates")
ax.legend()
fig.tight_layout()

OUTDIR.mkdir(exist_ok=True)
X.to_csv(OUTDIR / "winn_normalized_data.csv")
cv_control_data.to_csv(OUTDIR / "cv_control_data.csv", index=False)
cv_all_data.to_csv(OUTDIR / "cv_all_data.csv", index=False)
cv_combined_data.to_csv(OUTDIR / "cv_combined_data.csv", index=False)
metabolite_intensities.to_csv(OUTDIR / "original_metabolite_data.csv")
sample_metadata.to_csv(OUTDIR / "sample_metadata.csv", index=False)
fig.savefig(OUTDIR / "cv_plot.png", dpi=300)

print(f"✅ Done. Files saved in: {OUTDIR.resolve()}")
