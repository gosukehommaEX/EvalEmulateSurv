# =============================================================================
# test-EvalEmulateSurv.R
#
# Tests for the EvalEmulateSurv package.
# Run with: devtools::test()
# =============================================================================

library(EvalEmulateSurv)

# -----------------------------------------------------------------------------
# Shared test fixtures
# -----------------------------------------------------------------------------

# Full dummy IPD (both endpoints, both groups)
ipd <- gen_dummy_data(seed = 42)

# OS endpoint
ipd_os  <- ipd[ipd$TYPE == "OS", ]
grp1_os <- ipd_os[ipd_os$GROUP == "Drug",    ]
grp2_os <- ipd_os[ipd_os$GROUP == "Control", ]

# PFS endpoint
ipd_pfs  <- ipd[ipd$TYPE == "PFS", ]
grp1_pfs <- ipd_pfs[ipd_pfs$GROUP == "Drug",    ]
grp2_pfs <- ipd_pfs[ipd_pfs$GROUP == "Control", ]

tau_os  <- 36
tau_pfs <- 24


# =============================================================================
# gen_dummy_data()
# =============================================================================

test_that("gen_dummy_data: returns correct structure", {
  d <- gen_dummy_data(seed = 1)
  expect_s3_class(d, "data.frame")
  expect_named(d, c("ID", "GROUP", "TYPE", "TIME", "EVENT", "SEX", "REGION"))
})

test_that("gen_dummy_data: correct number of rows", {
  d <- gen_dummy_data(n_drug = 50, n_control = 60, seed = 1)
  # (50 + 60) * 2 endpoints = 220 rows
  expect_equal(nrow(d), 220L)
})

test_that("gen_dummy_data: GROUP levels are Drug and Control", {
  d <- gen_dummy_data(seed = 1)
  expect_setequal(unique(d$GROUP), c("Drug", "Control"))
})

test_that("gen_dummy_data: TYPE levels are OS and PFS", {
  d <- gen_dummy_data(seed = 1)
  expect_setequal(unique(d$TYPE), c("OS", "PFS"))
})

test_that("gen_dummy_data: TIME is strictly positive", {
  d <- gen_dummy_data(seed = 1)
  expect_true(all(d$TIME > 0))
})

test_that("gen_dummy_data: EVENT is 0 or 1", {
  d <- gen_dummy_data(seed = 1)
  expect_true(all(d$EVENT %in% c(0L, 1L)))
})

test_that("gen_dummy_data: SEX levels are M and F", {
  d <- gen_dummy_data(seed = 1)
  expect_true(all(d$SEX %in% c("M", "F")))
})

test_that("gen_dummy_data: REGION levels are US, EU, OTHERS", {
  d <- gen_dummy_data(seed = 1)
  expect_true(all(d$REGION %in% c("US", "EU", "OTHERS")))
})

test_that("gen_dummy_data: reproducible with same seed", {
  d1 <- gen_dummy_data(seed = 99)
  d2 <- gen_dummy_data(seed = 99)
  expect_identical(d1, d2)
})

test_that("gen_dummy_data: covariates consistent within patient across endpoints", {
  d <- gen_dummy_data(seed = 42)
  # For each ID x GROUP, SEX and REGION should be identical across OS and PFS
  os_cov  <- d[d$TYPE == "OS",  c("ID", "GROUP", "SEX", "REGION")]
  pfs_cov <- d[d$TYPE == "PFS", c("ID", "GROUP", "SEX", "REGION")]
  merged  <- merge(os_cov, pfs_cov,
                   by = c("ID", "GROUP"), suffixes = c("_os", "_pfs"))
  expect_true(all(merged$SEX_os    == merged$SEX_pfs))
  expect_true(all(merged$REGION_os == merged$REGION_pfs))
})


# =============================================================================
# load_ipd()
# =============================================================================

test_that("load_ipd: accepts data.frame and returns data.frame", {
  out <- load_ipd(grp1_os)
  expect_s3_class(out, "data.frame")
  expect_false(inherits(out, "tbl_df"))
})

test_that("load_ipd: accepts tibble and returns data.frame", {
  skip_if_not_installed("dplyr")
  tbl <- dplyr::as_tibble(grp1_os)
  out <- load_ipd(tbl)
  expect_s3_class(out, "data.frame")
  expect_false(inherits(out, "tbl_df"))
})

test_that("load_ipd: accepts csv file path", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(grp1_os, tmp, row.names = FALSE)
  out <- load_ipd(tmp)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(grp1_os))
  unlink(tmp)
})

test_that("load_ipd: accepts single-object RData file", {
  tmp <- tempfile(fileext = ".RData")
  grp1_save <- grp1_os
  save(grp1_save, file = tmp)
  out <- load_ipd(tmp)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(grp1_os))
  unlink(tmp)
})

test_that("load_ipd: multi-object RData requires rdata_object", {
  tmp <- tempfile(fileext = ".RData")
  obj_a <- grp1_os
  obj_b <- grp2_os
  save(obj_a, obj_b, file = tmp)
  expect_error(load_ipd(tmp), regexp = "rdata_object")
  out <- load_ipd(tmp, rdata_object = "obj_a")
  expect_s3_class(out, "data.frame")
  unlink(tmp)
})

test_that("load_ipd: errors on non-existent file", {
  expect_error(load_ipd("no_such_file.csv"), regexp = "not found")
})

test_that("load_ipd: errors on unsupported extension", {
  tmp <- tempfile(fileext = ".xlsx")
  file.create(tmp)
  expect_error(load_ipd(tmp), regexp = "Unsupported")
  unlink(tmp)
})


# =============================================================================
# km_est()
# =============================================================================

test_that("km_est: returns correct list structure", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_named(res, c("km1", "km2", "tau", "t_max_obs", "time_var",
                      "event_var", "desc1", "desc2", "conf_level"))
  expect_s3_class(res$km1, "data.frame")
  expect_s3_class(res$km2, "data.frame")
  expect_named(res$km1, c("time", "surv", "var_surv", "surv_lo", "surv_hi"))
})

test_that("km_est: var_surv is non-negative", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(all(res$km1$var_surv >= 0))
  expect_true(all(res$km2$var_surv >= 0))
})

test_that("km_est: surv_lo <= surv <= surv_hi", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(all(res$km1$surv_lo <= res$km1$surv + 1e-10))
  expect_true(all(res$km1$surv    <= res$km1$surv_hi + 1e-10))
})

test_that("km_est: desc1 and desc2 have correct structure", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_named(res$desc1, c("n", "n_events", "n_censored", "event_rate",
                            "median", "ci_lower", "ci_upper"))
  expect_named(res$desc2, c("n", "n_events", "n_censored", "event_rate",
                            "median", "ci_lower", "ci_upper"))
})

test_that("km_est: desc1 counts are consistent", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_equal(res$desc1$n, nrow(grp1_os))
  expect_equal(res$desc1$n_events + res$desc1$n_censored, res$desc1$n)
  expect_equal(res$desc2$n, nrow(grp2_os))
})

test_that("km_est: event_rate is in [0, 1]", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(res$desc1$event_rate >= 0 & res$desc1$event_rate <= 1)
  expect_true(res$desc2$event_rate >= 0 & res$desc2$event_rate <= 1)
})

test_that("km_est: ci_lower <= median <= ci_upper", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  if (!is.na(res$desc1$median) && !is.na(res$desc1$ci_lower) &&
      !is.na(res$desc1$ci_upper)) {
    expect_true(res$desc1$ci_lower <= res$desc1$median + 1e-10)
    expect_true(res$desc1$median   <= res$desc1$ci_upper + 1e-10)
  }
})

test_that("km_est: agrees with survival::survfit within machine tolerance", {
  skip_if_not_installed("survival")
  d  <- grp1_os
  sf <- survival::survfit(
    survival::Surv(d$TIME, d$EVENT) ~ 1,
    data = d, conf.type = "log"
  )
  sf_times <- sf$time[sf$n.event > 0]
  sf_surv  <- sf$surv[sf$n.event > 0]
  sf_lo    <- sf$lower[sf$n.event > 0]
  sf_hi    <- sf$upper[sf$n.event > 0]

  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT",
                tau = max(d$TIME))
  km1 <- res$km1

  .eval <- function(km_df, col, t) {
    idx <- findInterval(t, km_df$time)
    if (idx == 0L) 1.0 else km_df[[col]][idx]
  }

  est_surv <- sapply(sf_times, .eval, km_df = km1, col = "surv")
  est_lo   <- sapply(sf_times, .eval, km_df = km1, col = "surv_lo")
  est_hi   <- sapply(sf_times, .eval, km_df = km1, col = "surv_hi")

  tol <- .Machine$double.eps^0.5
  expect_true(max(abs(est_surv - sf_surv)) < tol)
  expect_true(max(abs(est_lo   - sf_lo))   < tol)
  expect_true(max(abs(est_hi   - sf_hi))   < tol)
})

test_that("km_est: median agrees with survival::survfit within machine tolerance", {
  skip_if_not_installed("survival")
  d  <- grp1_os
  sf <- survival::survfit(
    survival::Surv(d$TIME, d$EVENT) ~ 1,
    data = d, conf.type = "log"
  )
  q   <- quantile(sf, probs = 0.5)
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT",
                tau = max(d$TIME))

  expect_equal(res$desc1$median,   unname(q$quantile), tolerance = 1e-10)
  expect_equal(res$desc1$ci_lower, unname(q$lower),    tolerance = 1e-10)
  expect_equal(res$desc1$ci_upper, unname(q$upper),    tolerance = 1e-10)
})

test_that("km_est: S(0) = 1 for both groups", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_equal(res$km1$surv[1], 1.0)
  expect_equal(res$km2$surv[1], 1.0)
})

test_that("km_est: survival is monotonically non-increasing", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(all(diff(res$km1$surv) <= 0))
  expect_true(all(diff(res$km2$surv) <= 0))
})

test_that("km_est: survival is in [0, 1]", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(all(res$km1$surv >= 0 & res$km1$surv <= 1))
  expect_true(all(res$km2$surv >= 0 & res$km2$surv <= 1))
})

test_that("km_est: tau stored correctly", {
  res <- km_est(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_equal(res$tau, tau_os)
})

test_that("km_est: agrees with survival::survfit within machine tolerance", {
  skip_if_not_installed("survival")
  # Build KM using survival::survfit
  d <- grp1_os
  sf <- survival::survfit(
    survival::Surv(d$TIME, d$EVENT) ~ 1,
    data = d
  )
  # Extract unique event times and corresponding survival from survfit
  sf_times <- sf$time[sf$n.event > 0]
  sf_surv  <- sf$surv[sf$n.event > 0]

  # km_est result evaluated at the same event times
  res     <- km_est(grp1_os, grp2_os,
                    time_var = "TIME", event_var = "EVENT",
                    tau = max(d$TIME))
  km1     <- res$km1
  eval_surv <- sapply(sf_times, function(tt) {
    idx <- findInterval(tt, km1$time)
    if (idx == 0L) 1.0 else km1$surv[idx]
  })

  expect_true(
    all(abs(eval_surv - sf_surv) < .Machine$double.eps^0.5),
    info = sprintf(
      "Max deviation from survfit: %.2e (tolerance: %.2e)",
      max(abs(eval_surv - sf_surv)),
      .Machine$double.eps^0.5
    )
  )
})

test_that("km_est: warns when tau exceeds t_max_obs in both groups", {
  expect_warning(
    km_est(grp1_os, grp2_os,
           time_var = "TIME", event_var = "EVENT",
           tau = 1e6),
    regexp = "tau"
  )
})


# =============================================================================
# nauc()
# =============================================================================

test_that("nauc: returns correct list structure (overall)", {
  res <- nauc(grp1_os, grp2_os,
              time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_named(res, c("by_stratum", "aggregate", "method_avg", "tau"))
  expect_s3_class(res$by_stratum, "data.frame")
  expect_s3_class(res$aggregate,  "data.frame")
})

test_that("nauc: value is in [0, 1]", {
  res <- nauc(grp1_os, grp2_os,
              time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(res$by_stratum$nauc >= 0 & res$by_stratum$nauc <= 1)
})

test_that("nauc: identical data yields 0", {
  res <- nauc(grp1_os, grp1_os,
              time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_equal(res$by_stratum$nauc, 0, tolerance = 1e-10)
})

test_that("nauc: subgroup results have correct subgroup_spec column", {
  res <- nauc(grp1_os, grp2_os,
              time_var     = "TIME", event_var = "EVENT",
              tau          = tau_os,
              subgroup_var = "REGION")
  expect_true(all(res$by_stratum$subgroup_spec == "REGION"))
  expect_setequal(res$by_stratum$stratum,
                  sort(unique(c(grp1_os$REGION, grp2_os$REGION))))
})

test_that("nauc: interaction subgroup_var produces correct number of strata", {
  res <- nauc(grp1_os, grp2_os,
              time_var     = "TIME", event_var = "EVENT",
              tau          = tau_os,
              subgroup_var = "SEX:REGION")
  n_expected <- length(unique(
    interaction(grp1_os[, c("SEX", "REGION")], sep = " / ")
  ))
  expect_equal(nrow(res$by_stratum), n_expected)
})

test_that("nauc: include_overall prepends Overall row", {
  res <- nauc(grp1_os, grp2_os,
              time_var        = "TIME", event_var = "EVENT",
              tau             = tau_os,
              subgroup_var    = "REGION",
              include_overall = TRUE)
  expect_true("Overall" %in% res$by_stratum$subgroup_spec)
  expect_equal(res$by_stratum$subgroup_spec[1], "Overall")
})

test_that("nauc: aggregate equals overall nauc when subgroup_var = NULL", {
  res <- nauc(grp1_os, grp2_os,
              time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_equal(res$aggregate$nauc_aggregate, res$by_stratum$nauc)
})


# =============================================================================
# ksdist()
# =============================================================================

test_that("ksdist: returns correct list structure (overall)", {
  res <- ksdist(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_named(res, c("by_stratum", "aggregate", "method_avg", "tau"))
  expect_named(res$by_stratum,
               c("subgroup_spec", "stratum", "n_group1", "n_group2",
                 "n_total", "ks", "t_max_gap"))
})

test_that("ksdist: value is in [0, 1]", {
  res <- ksdist(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(res$by_stratum$ks >= 0 & res$by_stratum$ks <= 1)
})

test_that("ksdist: identical data yields 0", {
  res <- ksdist(grp1_os, grp1_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_equal(res$by_stratum$ks, 0, tolerance = 1e-10)
})

test_that("ksdist: t_max_gap is within [0, tau]", {
  res <- ksdist(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_true(res$by_stratum$t_max_gap >= 0)
  expect_true(res$by_stratum$t_max_gap <= tau_os)
})

test_that("ksdist: subgroup results have correct subgroup_spec column", {
  res <- ksdist(grp1_os, grp2_os,
                time_var     = "TIME", event_var = "EVENT",
                tau          = tau_os,
                subgroup_var = "SEX")
  expect_true(all(res$by_stratum$subgroup_spec == "SEX"))
  expect_setequal(res$by_stratum$stratum,
                  sort(unique(c(grp1_os$SEX, grp2_os$SEX))))
})

test_that("ksdist: include_overall prepends Overall row", {
  res <- ksdist(grp1_os, grp2_os,
                time_var        = "TIME", event_var = "EVENT",
                tau             = tau_os,
                subgroup_var    = "SEX",
                include_overall = TRUE)
  expect_equal(res$by_stratum$subgroup_spec[1], "Overall")
})

test_that("ksdist: KS <= 1 for interaction subgroup", {
  res <- ksdist(grp1_os, grp2_os,
                time_var     = "TIME", event_var = "EVENT",
                tau          = tau_os,
                subgroup_var = "SEX:REGION")
  expect_true(all(res$by_stratum$ks >= 0 & res$by_stratum$ks <= 1))
})


# =============================================================================
# rmse()
# =============================================================================

test_that("rmse: returns correct list structure (overall)", {
  res <- rmse(grp1_os, grp2_os,
              time_var = "TIME", event_var = "EVENT")
  expect_named(res, c("by_stratum", "aggregate", "method_avg",
                      "tau_input", "n_points"))
  expect_equal(res$n_points, 100L)
})

test_that("rmse: value is non-negative", {
  res <- rmse(grp1_os, grp2_os,
              time_var = "TIME", event_var = "EVENT")
  expect_true(res$by_stratum$rmse >= 0)
})

test_that("rmse: identical data yields 0", {
  res <- rmse(grp1_os, grp1_os,
              time_var = "TIME", event_var = "EVENT")
  expect_equal(res$by_stratum$rmse, 0, tolerance = 1e-10)
})

test_that("rmse: tau_input stored correctly", {
  res_null <- rmse(grp1_os, grp2_os,
                   time_var = "TIME", event_var = "EVENT", tau = NULL)
  res_tau  <- rmse(grp1_os, grp2_os,
                   time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_null(res_null$tau_input)
  expect_equal(res_tau$tau_input, tau_os)
})

test_that("rmse: n_points argument is respected", {
  res <- rmse(grp1_os, grp2_os,
              time_var = "TIME", event_var = "EVENT", n_points = 50L)
  expect_equal(res$n_points, 50L)
})

test_that("rmse: subgroup results all non-negative", {
  res <- rmse(grp1_os, grp2_os,
              time_var     = "TIME", event_var = "EVENT",
              subgroup_var = "REGION")
  expect_true(all(res$by_stratum$rmse >= 0))
})

test_that("rmse: include_overall prepends Overall row", {
  res <- rmse(grp1_os, grp2_os,
              time_var        = "TIME", event_var = "EVENT",
              subgroup_var    = "SEX",
              include_overall = TRUE)
  expect_equal(res$by_stratum$subgroup_spec[1], "Overall")
})


# =============================================================================
# kmplot()
# =============================================================================

test_that("kmplot: returns a named list", {
  res <- kmplot(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_type(res, "list")
  expect_named(res, "Overall")
})

test_that("kmplot: each element is a ggplot object", {
  res <- kmplot(grp1_os, grp2_os,
                time_var = "TIME", event_var = "EVENT", tau = tau_os)
  expect_s3_class(res[["Overall"]], "gg")
})

test_that("kmplot: subgroup_var produces correctly named list element", {
  res <- kmplot(grp1_os, grp2_os,
                time_var     = "TIME", event_var = "EVENT",
                tau          = tau_os,
                subgroup_var = "REGION",
                ncol         = 3)
  expect_true("REGION" %in% names(res))
  expect_s3_class(res[["REGION"]], "gg")
})

test_that("kmplot: interaction subgroup_var produces ggplot", {
  res <- kmplot(grp1_os, grp2_os,
                time_var     = "TIME", event_var = "EVENT",
                tau          = tau_os,
                subgroup_var = "SEX:REGION",
                ncol         = 3)
  expect_s3_class(res[["SEX:REGION"]], "gg")
})

test_that("kmplot: include_overall adds Overall panel to subgroup plot", {
  res <- kmplot(grp1_os, grp2_os,
                time_var        = "TIME", event_var = "EVENT",
                tau             = tau_os,
                subgroup_var    = "SEX",
                include_overall = TRUE,
                ncol            = 2)
  p <- res[["SEX"]]
  # Check that "Overall" appears in the facet labels
  built   <- ggplot2::ggplot_build(p)
  facet_labels <- unlist(built$layout$layout[["stratum"]])
  expect_true("Overall" %in% facet_labels)
})

test_that("kmplot: multiple subgroup_var returns multiple list elements", {
  res <- kmplot(grp1_os, grp2_os,
                time_var     = "TIME", event_var = "EVENT",
                tau          = tau_os,
                subgroup_var = c("SEX", "REGION"),
                ncol         = 2)
  expect_named(res, c("SEX", "REGION"))
})


# =============================================================================
# summarytable()
# =============================================================================

test_that("summarytable: returns EvalMetrics class", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  expect_s3_class(res, "EvalMetrics")
  expect_s3_class(res, "data.frame")
})

test_that("summarytable: required columns are present", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  expect_true(all(c("subgroup_spec", "stratum",
                    "n_group1", "n_group2",
                    "events_group1", "events_group2",
                    "median_group1", "median_ci_group1",
                    "median_group2", "median_ci_group2",
                    "nauc", "ks", "t_max_gap", "rmse",
                    "nauc_aggregate", "ks_aggregate",
                    "rmse_aggregate") %in% names(res)))
})

test_that("summarytable: event counts are consistent with n", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  ov  <- res[res$subgroup_spec == "Overall", ]
  expect_true(ov$events_group1 <= ov$n_group1)
  expect_true(ov$events_group2 <= ov$n_group2)
})

test_that("summarytable: median_ci columns are character strings", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  expect_type(res$median_ci_group1, "character")
  expect_type(res$median_ci_group2, "character")
})

test_that("summarytable: label_group1/2 stored in attributes", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var     = "TIME", event_var = "EVENT",
                      tau          = tau_os,
                      label_group1 = "Original",
                      label_group2 = "Emulated")
  expect_equal(attr(res, "label_group1"), "Original")
  expect_equal(attr(res, "label_group2"), "Emulated")
})

test_that("summarytable: Overall block is always present", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  expect_true("Overall" %in% res$subgroup_spec)
})

test_that("summarytable: subgroup block present when specified", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var      = "TIME", event_var = "EVENT",
                      tau           = tau_os,
                      subgroup_vars = "REGION")
  expect_true("REGION" %in% res$subgroup_spec)
})

test_that("summarytable: interaction subgroup_var produces correct spec", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var      = "TIME", event_var = "EVENT",
                      tau           = tau_os,
                      subgroup_vars = "SEX:REGION")
  expect_true("SEX:REGION" %in% res$subgroup_spec)
})

test_that("summarytable: nauc_aggregate is NA for Overall block", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  ov_rows <- res[res$subgroup_spec == "Overall", ]
  expect_true(all(is.na(ov_rows$nauc_aggregate)))
  expect_true(all(is.na(ov_rows$ks_aggregate)))
  expect_true(all(is.na(ov_rows$rmse_aggregate)))
})

test_that("summarytable: all metric values are non-negative", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var      = "TIME", event_var = "EVENT",
                      tau           = tau_os,
                      subgroup_vars = c("SEX", "REGION"))
  expect_true(all(res$nauc >= 0, na.rm = TRUE))
  expect_true(all(res$ks   >= 0, na.rm = TRUE))
  expect_true(all(res$rmse >= 0, na.rm = TRUE))
})

test_that("summarytable: nauc and ks are in [0, 1]", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var      = "TIME", event_var = "EVENT",
                      tau           = tau_os,
                      subgroup_vars = "REGION")
  expect_true(all(res$nauc <= 1, na.rm = TRUE))
  expect_true(all(res$ks   <= 1, na.rm = TRUE))
})

test_that("summarytable: metadata attributes are attached", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  expect_equal(attr(res, "tau"),        tau_os)
  expect_equal(attr(res, "time_var"),   "TIME")
  expect_equal(attr(res, "event_var"),  "EVENT")
  expect_equal(attr(res, "n_points"),   100L)
})

test_that("summarytable: print method runs without error (plain text)", {
  res <- summarytable(grp1_os, grp2_os,
                      time_var      = "TIME", event_var = "EVENT",
                      tau           = tau_os,
                      subgroup_vars = "REGION")
  expect_output(print(res), regexp = "EvalEmulateSurv")
})

test_that("summarytable: print method with kable = TRUE runs without error", {
  skip_if_not_installed("knitr")
  res <- summarytable(grp1_os, grp2_os,
                      time_var = "TIME", event_var = "EVENT",
                      tau = tau_os)
  expect_no_error(print(res, kable = TRUE, kable_format = "pipe"))
})
