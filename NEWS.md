# winn 0.1.2

- Added `winn_ablation()` for fixed-parameter cumulative and gate-factorial
  analyses using the same correction operators as `winn()`.
- Added optional decision and correction diagnostics to the drift, batch, and
  dilution-normalization functions. Their default matrix-returning behavior is
  unchanged.
- Made selective, forced-all, and disabled drift and batch gates explicit for
  controlled method evaluation.
- Documented the directional, upper-priority single-tail MAD policy and added
  regression tests for upper-only, lower-only, simultaneous-tail, and
  zero-MAD inputs.
- Added tests for ablation equivalence, disabled-stage pass-through behavior,
  diagnostic structure, and standard-smoother fallback reporting.
