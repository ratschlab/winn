import argparse, json
import pandas as pd
from .winn import winn

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xlsx", required=True)
    p.add_argument("--sheet", default="intensities1")
    p.add_argument("--meta", required=True)
    p.add_argument("--test", default="Ljung-Box", choices=["Ljung-Box","DW"])
    p.add_argument("--detrend", default="mean", choices=["mean","spline"])
    p.add_argument("--spline-method", default="conservative", choices=["conservative","standard"])
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--median", default="normalize", choices=["normalize","shrink","none"])
    p.add_argument("--batch-method", default="combat", choices=["combat","anova"])
    p.add_argument("--lag", type=int, default=20)
    p.add_argument("--scale-by-batch", action="store_true")
    p.add_argument("--params", default="fixed", choices=["fixed","auto"])
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
        parameters=args.params,
        fdr_threshold=args.fdr,
        median_adjustment=args.median,
        detrend_non_autocorrelated=args.detrend,
        remove_batch_effects=args.batch_method,
        test=args.test,
        lag=args.lag,
        scale_by_batch=args.scale_by_batch,
        spline_method=args.spline_method,
    )
    out.to_csv("winn_normalized_data.csv")
    print("Saved: winn_normalized_data.csv")
