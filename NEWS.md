# winn 0.1.4

- Made drift testing and smoothing honor `run_order` independently of matrix
  column order. Values are sorted within batch for testing and correction and
  restored to the caller's original column order.
- Changed the Ljung-Box default to `fitdf = 0`, because WiNN does not first fit
  an ARMA model. The legacy `fitdf = 1` behavior remains available explicitly
  through `ljung_box_fitdf = 1`.
- Added strict validation for finite, non-negative log-scale inputs; complete
  batch and run-order metadata; unique within-batch run orders; control
  indices; dilution factors; and within-batch scaling.
- Documented the preprocessing policy for zero and censored observations.
  Numeric zeros remain observed values and are never silently converted to
  missing values.
- Added `return_details = TRUE` to `winn()` for machine-readable selected
  parameters, batch assignments, and per-stage decisions while retaining the
  corrected-matrix default.
- Clarified fixed versus auto-mode settings and the selective ANOVA versus
  global ComBat correction scopes.
- Removed interactive test code from package sources and added release CI.

# winn 0.1.3

- Replaced the legacy upper-priority MAD branch with independent two-tail
  shrinkage. Each metabolite's median, MAD, thresholds, and upper/lower masks
  are computed once from the original finite values; both tails are then
  adjusted in one non-iterative pass.
- Added regression coverage for simultaneous-tail adjustment, one-pass
  threshold reuse, missing/non-finite values, and zero-MAD features.

# winn 0.1.2

- Added `winn_ablation()` for fixed-parameter cumulative and gate-factorial
  analyses using the same correction operators as `winn()`.
- Added optional decision and correction diagnostics to the drift, batch, and
  dilution-normalization functions. Their default matrix-returning behavior is
  unchanged.
- Made selective, forced-all, and disabled drift and batch gates explicit for
  controlled method evaluation.
- Documented the historical directional MAD behavior and added regression
  tests that made the legacy branch explicit before its replacement in 0.1.3.
- Added tests for ablation equivalence, disabled-stage pass-through behavior,
  diagnostic structure, and standard-smoother fallback reporting.
