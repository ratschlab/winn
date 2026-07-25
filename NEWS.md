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
