# EvalEmulateSurv 0.1.0

## Initial release

### New functions

* `gen_dummy_data()` -- generates dummy individual patient data (IPD) for
  examples and testing, with survival times drawn from exponential
  distributions and covariates SEX and REGION.

* `load_ipd()` -- loads IPD from a `data.frame`, tibble, `.csv` file, or
  single-object `.RData` / `.rda` file into a plain `data.frame`.

* `km_est()` -- computes Kaplan-Meier step functions for two groups without
  relying on the `survival` package.

* `nauc()` -- computes the Normalised Area Under the Curve (NAUC) between
  two KM estimates, with support for subgroup variables and N-way interaction
  terms.

* `ksdist()` -- computes the Kolmogorov-Smirnov (K-S) distance between two
  KM estimates, with support for subgroup variables and N-way interaction
  terms.

* `rmse()` -- computes the Root Mean Square Error (RMSE) between two KM
  estimates evaluated at evenly spaced time points, with support for subgroup
  variables and N-way interaction terms.

* `kmplot()` -- draws Kaplan-Meier plots with NAUC, K-S distance, and RMSE
  annotations using `facet_wrap()`, supporting subgroup variables and N-way
  interaction terms.

* `summarytable()` -- computes all three metrics for the overall population
  and subgroup variables in a single call. Returns an `EvalMetrics` S3 object
  with a formatted `print()` method supporting plain text, `kable`, and HTML
  output via `kableExtra`.
