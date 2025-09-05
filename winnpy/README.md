
# winnpy — Quick Usage

**WINN for Python**: end-to-end correction for metabolomics matrices
(rows = metabolites, columns = samples).

Pipeline: MAD outliers → drift correction (by run order) → batch correction (ANOVA/ComBat) → median/PQN → optional per-batch scaling. Processing is done in **log1p** space and returned with `expm1`.

---

## 1) Install (from this repo)

> Python 3.10+

```powershell
# Windows PowerShell (run at repo root)
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -U pip
pip install -e ".[ruptures,combat]"   # optional extras; omit if not needed
```

```bash
# macOS / Linux
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e ".[ruptures,combat]"
```

Sanity check:

```bash
python - <<'PY'
from winnpy import winn
print("winn imported OK")
PY
```

---

## 2) One-command run (provided driver)

Put files in `data/`:

* `data/raetsch1_DATA_FILTERED.xlsx` (sheet `intensities1`)
* `data/per_sample_info.csv`

Run:

```bash
python scripts/run_winn_package.py
```

Outputs (in `scripts/results/`):

* `winn_normalized_data.csv`
* `cv_control_data.csv`, `cv_all_data.csv`, `cv_combined_data.csv`
* `original_metabolite_data.csv`, `sample_metadata.csv`
* `cv_plot.png`


---

## 3) Common gotchas

* **Module not found** → activate the env where you installed:

  * Windows: `.\.venv\Scripts\Activate.ps1`
  * macOS/Linux: `source .venv/bin/activate`
* **ComBat not available** → set `remove_batch_effects="anova"` or `pip install combat`.

---
