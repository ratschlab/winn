from pathlib import Path
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

# === 1) 路径设置 ===
ROOT = Path(__file__).resolve().parent        # scripts/
DATA_DIR = ROOT.parent / "data"               # Pywinn/data

XLSX_PATH = DATA_DIR / "raetsch1_DATA_FILTERED.xlsx"
META_PATH = DATA_DIR / "per_sample_info.csv"

OUTDIR = ROOT / "results"                     # scripts/results
OUTDIR.mkdir(exist_ok=True)


# === 2) 读数据（保持和原脚本一致） ===
metabolite_intensities = pd.read_excel(XLSX_PATH, sheet_name="intensities1")
if not np.issubdtype(metabolite_intensities.dtypes[0], np.number):
    metabolite_intensities = metabolite_intensities.set_index(metabolite_intensities.columns[0])
metabolite_intensities = metabolite_intensities.apply(pd.to_numeric, errors="coerce")

sample_metadata = pd.read_csv(META_PATH)

# 修正列名（和你原脚本一致）
metabolite_intensities.columns = (
    metabolite_intensities.columns.astype(str)
      .str.replace("HG20.535", "HG2020.535", regex=False)
      .str.replace("HG20.709", "HG2020.709", regex=False)
)

# 按 sampleID 重排
if "sampleID" not in sample_metadata.columns:
    raise ValueError("Metadata is missing 'sampleID'")
cols = sample_metadata["sampleID"].astype(str).tolist()
metabolite_intensities.columns = metabolite_intensities.columns.astype(str)
missing = set(cols) - set(metabolite_intensities.columns)
if missing:
    raise ValueError(f"These samples are not found: {sorted(missing)}")
metabolite_intensities = metabolite_intensities.loc[:, cols]

# === 3) 组装 WiNN 参数并调用你打包的 Python 包 ===
# 依赖：pip install -e .[combat,ruptures]（combat/ruptures 可选）
try:
    # 假设 __init__.py 已导出 winn；否则用 from winnpy.winn import winn
    from winnpy import winn
except Exception:
    from winnpy.winn import winn

# batch = as.numeric(as.factor(plate))；run_order = order；control = sample=="control"
if "plate" not in sample_metadata.columns:
    raise ValueError("Metadata is missing 'plate'")
if "order" not in sample_metadata.columns:
    raise ValueError("Metadata is missing 'order'")
if "sample" not in sample_metadata.columns:
    raise ValueError("Metadata is missing 'sample' (used to identify controls)")

batch_codes = pd.Categorical(sample_metadata["plate"]).codes + 1
run_order = sample_metadata["order"].to_numpy()
control_cols = sample_metadata.loc[
    sample_metadata["sample"].astype(str).str.lower() == "control", "sampleID"
].astype(str).tolist() or None

# 与你们 R 流程一致的固定参数（若没装 combat，会自动降级为 anova）
params = dict(
    batch=batch_codes,
    run_order=run_order,
    control_samples=control_cols,
    parameters="fixed",
    fdr_threshold=1.0,
    median_adjustment="normalize",
    detrend_non_autocorrelated="mean",
    remove_batch_effects="combat",  # 尝试 ComBat
    test="Ljung-Box",
    lag=20,
    scale_by_batch=False,
)

try:
    X = winn(data=metabolite_intensities, **params)
except ImportError as e:
    # 没装 combat：降级为 anova
    params["remove_batch_effects"] = "anova"
    X = winn(data=metabolite_intensities, **params)

# === 4) 计算 CV（control / all），与原脚本一致 ===
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

# === 5) 画图（均值 ± SE） ===
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

# === 6) 落盘（文件名保持和你原脚本一致） ===
X.to_csv(OUTDIR / "winn_normalized_data.csv")
cv_control_data.to_csv(OUTDIR / "cv_control_data.csv", index=False)
cv_all_data.to_csv(OUTDIR / "cv_all_data.csv", index=False)
cv_combined_data.to_csv(OUTDIR / "cv_combined_data.csv", index=False)
metabolite_intensities.to_csv(OUTDIR / "original_metabolite_data.csv")
sample_metadata.to_csv(OUTDIR / "sample_metadata.csv", index=False)
fig.savefig(OUTDIR / "cv_plot.png", dpi=300)

print(f"✅ Done. Files saved in: {OUTDIR.resolve()}")
