#' Kaplan-Meier Estimation for Two Groups
#'
#' @description
#' Computes Kaplan-Meier (KM) step functions for two groups from raw
#' individual patient data (IPD). The KM estimator is implemented from
#' scratch without relying on the \pkg{survival} package. The output is
#' a structured list that serves as the common input for \code{nauc()},
#' \code{ksdist()}, \code{rmse()}, \code{kmplot()}, and
#' \code{summarytable()}. Input data can be supplied as a \code{data.frame},
#' tibble, path to a \code{.csv} file, or path to an \code{.RData} /
#' \code{.rda} file containing a single object, via \code{\link{load_ipd}}.
#' For \code{.RData} files containing multiple objects, call
#' \code{\link{load_ipd}} explicitly with the \code{rdata_object} argument
#' before passing the result here.
#'
#' @details
#' The KM estimate at each event time \eqn{t_{(j)}} is computed as:
#' \deqn{
#'   \hat{S}(t_{(j)}) = \prod_{i:\, t_{(i)} \le t_{(j)}}
#'   \left(1 - \frac{d_i}{n_i}\right)
#' }
#' where \eqn{d_i} is the number of events and \eqn{n_i} is the number
#' at risk at time \eqn{t_{(i)}}. The number at risk is obtained via
#' \code{findInterval()} on the sorted observation times, avoiding any
#' dependency on external survival packages.
#'
#' The returned step function represents \eqn{\hat{S}(t)} on the
#' half-open interval \eqn{[t_{(j)}, t_{(j+1)})}, i.e., the survival
#' probability is constant between consecutive event times and drops
#' at each event time (left-continuous convention).
#'
#' \strong{Greenwood variance and pointwise confidence band:}
#'
#' The variance of \eqn{\hat{S}(t)} is estimated by Greenwood's formula:
#' \deqn{
#'   \widehat{\mathrm{Var}}(\hat{S}(t)) =
#'   \hat{S}(t)^2 \sum_{t_{(j)} \le t}
#'   \frac{d_j}{n_j (n_j - d_j)}
#' }
#' A pointwise \code{conf_level} confidence band is constructed on the log
#' scale (same as \code{survival::survfit} with \code{conf.type = "log"}):
#' \deqn{
#'   \exp\!\left(\log \hat{S}(t) \pm
#'   z_{\alpha/2} \frac{\sqrt{\widehat{\mathrm{Var}}(\hat{S}(t))}}
#'   {\hat{S}(t)}\right)
#' }
#'
#' \strong{Median survival time and CI:}
#'
#' The median is the smallest \eqn{t} such that \eqn{\hat{S}(t) \le 0.5}.
#' The CI for the median is the range of times where the pointwise confidence
#' band crosses 0.5: \code{ci_lower} is the smallest \eqn{t} where the lower
#' band \eqn{\le 0.5}, and \code{ci_upper} is the smallest \eqn{t} where the
#' upper band \eqn{\le 0.5}. \code{NA} is returned when a bound is not
#' estimable within the observed range.
#'
#' @param data_group1 Input for group 1 (e.g., original IPD). Accepted
#'   formats: a \code{data.frame} or tibble already in memory, a character
#'   path to a \code{.csv} file, or a character path to an \code{.RData} /
#'   \code{.rda} file containing a single object. Must contain columns
#'   specified by \code{time_var} and \code{event_var}.
#'   See \code{\link{load_ipd}} for details.
#' @param data_group2 Input for group 2 (e.g., emulated IPD). Same accepted
#'   formats as \code{data_group1}. Must contain columns specified by
#'   \code{time_var} and \code{event_var}.
#' @param time_var  Character. Name of the time-to-event column.
#'   Default \code{"TIME"}.
#' @param event_var Character. Name of the event indicator column
#'   (1 = event, 0 = censored). Default \code{"EVENT"}.
#' @param tau Numeric scalar. Truncation time. KM estimates are
#'   computed up to \code{tau}. If \code{NULL} (default), the maximum
#'   observed time across both groups is used.
#' @param conf_level Numeric in \code{(0, 1)}. Confidence level for the
#'   pointwise CI band and the median CI. Default \code{0.95}.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\code{km1}}{A \code{data.frame} with columns \code{time},
#'     \code{surv}, \code{var_surv} (Greenwood variance), \code{surv_lo}
#'     (lower pointwise CI band), and \code{surv_hi} (upper pointwise CI
#'     band) for group 1.}
#'   \item{\code{km2}}{Same structure as \code{km1} for group 2.}
#'   \item{\code{tau}}{The truncation time actually used.}
#'   \item{\code{t_max_obs}}{The maximum observed time across both
#'     groups (before truncation).}
#'   \item{\code{time_var}}{The name of the time variable used.}
#'   \item{\code{event_var}}{The name of the event variable used.}
#'   \item{\code{desc1}}{A named list of descriptive statistics for group 1:
#'     \code{n}, \code{n_events}, \code{n_censored}, \code{event_rate},
#'     \code{median}, \code{ci_lower}, \code{ci_upper}.}
#'   \item{\code{desc2}}{Same structure as \code{desc1} for group 2.}
#'   \item{\code{conf_level}}{The confidence level used.}
#' }
#'
#' @seealso \code{\link{load_ipd}}, \code{\link{nauc}},
#'   \code{\link{ksdist}}, \code{\link{rmse}}, \code{\link{kmplot}},
#'   \code{\link{summarytable}}
#'
#' @examples
#' # -------------------------------------------------------------------
#' # Prepare dummy data with gen_dummy_data()
#' # -------------------------------------------------------------------
#' ipd <- gen_dummy_data(seed = 42)
#'
#' # Filter to OS endpoint and split by group
#' ipd_os  <- ipd[ipd$TYPE == "OS", ]
#' grp1_os <- ipd_os[ipd_os$GROUP == "1",    ]
#' grp2_os <- ipd_os[ipd_os$GROUP == "2", ]
#'
#' # Filter to PFS endpoint and split by group
#' ipd_pfs  <- ipd[ipd$TYPE == "PFS", ]
#' grp1_pfs <- ipd_pfs[ipd_pfs$GROUP == "1",    ]
#' grp2_pfs <- ipd_pfs[ipd_pfs$GROUP == "2", ]
#'
#' # -------------------------------------------------------------------
#' # 1. OS, overall population -- data.frame input
#' # -------------------------------------------------------------------
#' km_os <- km_est(
#'   data_group1 = grp1_os,
#'   data_group2 = grp2_os,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' head(km_os$km1)   # time, surv, var_surv, surv_lo, surv_hi
#' km_os$desc1       # n, n_events, n_censored, event_rate, median, CI
#' km_os$desc2
#' km_os$tau
#' km_os$t_max_obs
#'
#' # -------------------------------------------------------------------
#' # 2. PFS, overall population
#' # -------------------------------------------------------------------
#' km_pfs <- km_est(
#'   data_group1 = grp1_pfs,
#'   data_group2 = grp2_pfs,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 24
#' )
#' km_pfs$desc1
#'
#' # -------------------------------------------------------------------
#' # 3. CSV file input via load_ipd()
#' # -------------------------------------------------------------------
#' tmp1 <- tempfile(fileext = ".csv")
#' tmp2 <- tempfile(fileext = ".csv")
#' write.csv(grp1_os, tmp1, row.names = FALSE)
#' write.csv(grp2_os, tmp2, row.names = FALSE)
#' km_csv <- km_est(
#'   data_group1 = tmp1,
#'   data_group2 = tmp2,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' head(km_csv$km1)
#' unlink(c(tmp1, tmp2))
#'
#' # -------------------------------------------------------------------
#' # 4. RData file (single object) input via load_ipd()
#' # -------------------------------------------------------------------
#' tmp_rda1 <- tempfile(fileext = ".RData")
#' tmp_rda2 <- tempfile(fileext = ".RData")
#' save(grp1_os, file = tmp_rda1)
#' save(grp2_os, file = tmp_rda2)
#' km_rda <- km_est(
#'   data_group1 = tmp_rda1,
#'   data_group2 = tmp_rda2,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' head(km_rda$km1)
#' unlink(c(tmp_rda1, tmp_rda2))
#'
#' # -------------------------------------------------------------------
#' # 5. RData file (multiple objects): use load_ipd() explicitly
#' # -------------------------------------------------------------------
#' tmp_multi <- tempfile(fileext = ".RData")
#' save(grp1_os, grp2_os, file = tmp_multi)
#' d1 <- load_ipd(tmp_multi, rdata_object = "grp1_os")
#' d2 <- load_ipd(tmp_multi, rdata_object = "grp2_os")
#' km_multi <- km_est(
#'   data_group1 = d1,
#'   data_group2 = d2,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' head(km_multi$km1)
#' unlink(tmp_multi)
#'
#' @export
km_est <- function(data_group1,
                   data_group2,
                   time_var   = "TIME",
                   event_var  = "EVENT",
                   tau        = NULL,
                   conf_level = 0.95) {

  # ---------------------------------------------------------------------------
  # Step 1: Load input data via load_ipd()
  # ---------------------------------------------------------------------------
  data_group1 <- load_ipd(data_group1)
  data_group2 <- load_ipd(data_group2)

  # ---------------------------------------------------------------------------
  # Step 2: Input validation
  # ---------------------------------------------------------------------------
  if (!time_var %in% names(data_group1))
    stop(sprintf("'%s' not found in data_group1.", time_var))
  if (!time_var %in% names(data_group2))
    stop(sprintf("'%s' not found in data_group2.", time_var))
  if (!event_var %in% names(data_group1))
    stop(sprintf("'%s' not found in data_group1.", event_var))
  if (!event_var %in% names(data_group2))
    stop(sprintf("'%s' not found in data_group2.", event_var))
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1)
    stop("'conf_level' must be a single number in (0, 1).")

  z <- qnorm(1 - (1 - conf_level) / 2)

  # ---------------------------------------------------------------------------
  # Internal helper: compute KM step function with Greenwood variance and
  # pointwise confidence band from raw vectors.
  #
  # Returns a data.frame with columns:
  #   time     : event times starting from 0
  #   surv     : KM estimate S(t) on [time[i], time[i+1])
  #   var_surv : Greenwood variance of S(t)
  #   surv_lo  : lower pointwise CI bound (log transformation)
  #   surv_hi  : upper pointwise CI bound (log transformation)
  #
  # Only event times up to tau_max are included. Observations beyond
  # tau_max are treated as censored at tau_max.
  # ---------------------------------------------------------------------------
  .km_stepfun <- function(time_vec, event_vec, tau_max) {
    event_times <- sort(unique(time_vec[event_vec == 1L]))
    event_times <- event_times[event_times <= tau_max]

    # If no events, survival stays at 1 throughout
    if (length(event_times) == 0L) {
      return(data.frame(time     = c(0, tau_max),
                        surv     = c(1, 1),
                        var_surv = c(0, 0),
                        surv_lo  = c(1, 1),
                        surv_hi  = c(1, 1)))
    }

    n_total     <- length(time_vec)
    time_sorted <- sort(time_vec)

    surv      <- 1.0
    greenwood <- 0.0
    s_vec     <- numeric(length(event_times))
    v_vec     <- numeric(length(event_times))

    for (i in seq_along(event_times)) {
      tt      <- event_times[i]
      # Number at risk just before tt (left-open interval)
      n_risk  <- n_total - findInterval(tt, time_sorted, left.open = TRUE)
      # Number of events at tt
      n_event <- sum(time_vec[event_vec == 1L] == tt)
      surv    <- surv * (1.0 - n_event / n_risk)

      # Greenwood increment: skip when n_risk == n_event (avoids division by 0)
      if (n_risk > n_event)
        greenwood <- greenwood + n_event / (n_risk * (n_risk - n_event))

      s_vec[i] <- surv
      v_vec[i] <- surv^2 * greenwood
    }

    # Pointwise CI on log scale (same as survival::survfit conf.type = "log")
    se_vec <- sqrt(v_vec)
    s_safe <- pmax(s_vec, .Machine$double.eps)   # avoid log(0)
    lo_vec <- pmin(pmax(exp(log(s_safe) - z * se_vec / s_safe), 0), 1)
    hi_vec <- pmin(pmax(exp(log(s_safe) + z * se_vec / s_safe), 0), 1)

    # Row i: S(t) is constant on [event_times[i], event_times[i+1])
    # Include t = 0 with S(0) = 1 as the starting row
    data.frame(
      time     = c(0,      event_times),
      surv     = c(1,      s_vec),
      var_surv = c(0,      v_vec),
      surv_lo  = c(1,      lo_vec),
      surv_hi  = c(1,      hi_vec)
    )
  }

  # ---------------------------------------------------------------------------
  # Internal helper: extract median survival time and CI from a KM data.frame.
  #
  # median  : smallest t where surv    <= 0.5
  # ci_lower: smallest t where surv_lo <= 0.5  (lower band crosses 0.5 earlier)
  # ci_upper: smallest t where surv_hi <= 0.5  (upper band crosses 0.5 later)
  # ---------------------------------------------------------------------------
  .extract_median_ci <- function(km_df) {
    .first_cross <- function(vals, times) {
      idx <- which(vals <= 0.5)
      if (length(idx) == 0L) NA_real_ else times[idx[1L]]
    }
    list(
      median   = .first_cross(km_df$surv,    km_df$time),
      ci_lower = .first_cross(km_df$surv_lo, km_df$time),
      ci_upper = .first_cross(km_df$surv_hi, km_df$time)
    )
  }

  # ---------------------------------------------------------------------------
  # Internal helper: compute descriptive statistics for one group.
  # Uses the full observation range (not truncated to tau) so that the
  # median and CI reflect the complete follow-up.
  # ---------------------------------------------------------------------------
  .desc_one <- function(time_vec, event_vec) {
    n          <- length(time_vec)
    n_events   <- sum(event_vec)
    n_censored <- n - n_events
    event_rate <- n_events / n

    # Build KM over full range for median estimation
    km_full <- .km_stepfun(time_vec, event_vec, max(time_vec))
    mci     <- .extract_median_ci(km_full)

    list(
      n          = n,
      n_events   = n_events,
      n_censored = n_censored,
      event_rate = event_rate,
      median     = mci$median,
      ci_lower   = mci$ci_lower,
      ci_upper   = mci$ci_upper
    )
  }

  # ---------------------------------------------------------------------------
  # Step 3: Extract vectors and determine tau
  # ---------------------------------------------------------------------------
  t1 <- data_group1[[time_var]]
  t2 <- data_group2[[time_var]]
  e1 <- as.integer(data_group1[[event_var]])
  e2 <- as.integer(data_group2[[event_var]])

  t_max_obs <- max(c(t1, t2), na.rm = TRUE)

  if (is.null(tau)) {
    tau <- t_max_obs
  } else {
    if (!is.numeric(tau) || length(tau) != 1L || tau <= 0)
      stop("'tau' must be a single positive number.")
    # Warn only when tau exceeds the maximum observed time in BOTH groups.
    # If tau exceeds only one group's maximum (e.g. when rmse() uses the
    # overall max as tau for a subgroup subset), no warning is needed
    # because the KM step function handles the extrapolation correctly by
    # carrying forward the last observed survival probability.
    t_max_g1 <- max(t1, na.rm = TRUE)
    t_max_g2 <- max(t2, na.rm = TRUE)
    if (tau > t_max_g1 && tau > t_max_g2) {
      warning(sprintf(
        "'tau' (%.4f) exceeds the maximum observed time in both groups (%.4f).",
        tau, t_max_obs
      ))
    }
  }

  # ---------------------------------------------------------------------------
  # Step 4: Compute KM step functions (with Greenwood variance) up to tau
  # ---------------------------------------------------------------------------
  km1 <- .km_stepfun(t1, e1, tau)
  km2 <- .km_stepfun(t2, e2, tau)

  # ---------------------------------------------------------------------------
  # Step 5: Compute descriptive statistics over the full observation range
  # ---------------------------------------------------------------------------
  desc1 <- .desc_one(t1, e1)
  desc2 <- .desc_one(t2, e2)

  # ---------------------------------------------------------------------------
  # Step 6: Return structured list
  # ---------------------------------------------------------------------------
  list(
    km1        = km1,
    km2        = km2,
    tau        = tau,
    t_max_obs  = t_max_obs,
    time_var   = time_var,
    event_var  = event_var,
    desc1      = desc1,
    desc2      = desc2,
    conf_level = conf_level
  )
}
