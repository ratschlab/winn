import argparse, json
import pandas as pd
from .winn import winn

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xlsx", required=True)
    p.add_argument("--sheet", default="intensities1")
    p.add_argument("--meta", required=True)
    args = p.parse_args()

    intens = pd.read_excel(args.xlsx, sheet_name=args.sheet)
    if not pd.api.types.is_numeric_dtype(intens.iloc[:,0]):
        intens = intens.set_index(intens.columns[0])

    meta = pd.read_csv(args.meta)
    batch = pd.Categorical(meta["plate"]).codes + 1
    run_order = meta["order"].to_numpy()
    ctrl = meta.loc[meta["sample"].astype(str).str.lower()=="control","sampleID"].astype(str).tolist() or None

    out = winn(
        data=intens.loc[:, meta["sampleID"].astype(str)],
        batch=batch,
        run_order=run_order,
        control_samples=ctrl,
        parameters="fixed",
        fdr_threshold=1.0,
        median_adjustment="normalize",
        detrend_non_autocorrelated="mean",
        remove_batch_effects="combat",  # Will throw error if combat package is not available, can change to anova
        test="Ljung-Box",
        lag=20,
        scale_by_batch=False,
    )
    out.to_csv("winn_normalized_data.csv")
    print("Saved: winn_normalized_data.csv")
