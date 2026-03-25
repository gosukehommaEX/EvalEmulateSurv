# EvalEmulateSurv 0.1.0

## Initial release

### New functions

* `gen_dummy_data()` -- generates dummy individual patient data (IPD) for
  examples and testing, with survival times drawn from exponential
  distributions and covariates SEX and REGION.

* `load_ipd()` -- loads IPD from a `data.frame`, tibble, `.csv` file, or
  single-object `.RData` / `.rda` file into a plain `data.frame`.

* `km_est()` -- computes Kaplan-Meier step functions for two groups without
  relying on the `survival` package. Returns KM survival estimates with
  Greenwood variance and pointwise confidence bands (`surv_lo`, `surv_hi`),
  plus descriptive statistics per group (`desc1`, `desc2`): sample size,
  event count, censored count, event rate, median survival time, and
  its confidence interval (Greenwood + log transformation, equivalent to
  `survival::survfit` with `conf.type = "log"`).

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

* `summarytable()` -- computes all three metrics together with descriptive
  survival statistics (sample size, event count, median survival time with
  confidence interval) for the overall population and subgroup variables in a
  single call. Returns an `EvalMetrics` S3 object with a formatted `print()`
  method supporting plain text, `kable`, and HTML output via `kableExtra`.
