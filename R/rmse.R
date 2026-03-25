#' Root Mean Square Error (RMSE) Between Two KM Estimates
#'
#' @description
#' Computes the Root Mean Square Error (RMSE) between two Kaplan-Meier
#' survival curves as defined in Lang et al. (2025, RESOLVE-IPD) and
#' Zhang et al. (2024, SurvdigitizeR). The RMSE measures the average
#' discrepancy between the estimated survival probabilities of group 1
#' (ground truth) and group 2 (emulated data) at 100 evenly spaced time
#' points from 0 to the maximum observation time. KM estimation is
#' performed internally via \code{\link{km_est}}. Input data can be
#' supplied as a \code{data.frame}, tibble, path to a \code{.csv} file,
#' or path to an \code{.RData} / \code{.rda} file containing a single
#' object, via \code{\link{load_ipd}}. For \code{.RData} files containing
#' multiple objects, call \code{\link{load_ipd}} explicitly with the
#' \code{rdata_object} argument before passing the result here.
#'
#' @details
#' When \code{tau = NULL} (default), the upper limit of the evaluation grid
#' is determined per stratum from the subsetted data: each stratum uses
#' its own maximum observed time. This ensures that the per-stratum RMSE
#' returned here is consistent with the value displayed in the corresponding
#' \code{kmplot()} panel, which receives only the stratum-subset data.
#' When \code{tau} is explicitly specified, the same value is used for all
#' strata.
#'
#' Let \eqn{t_1, \ldots, t_{100}} be 100 evenly spaced time points from
#' \eqn{0} to \eqn{\tau}. The RMSE for a single stratum is:
#'
#' \deqn{
#'   \mathrm{RMSE} =
#'   \sqrt{\frac{1}{100} \sum_{k=1}^{100}
#'   \left( \hat{S}(t_k) - \tilde{S}(t_k) \right)^2}
#' }
#'
#' where \eqn{\hat{S}(t_k)} is the KM estimate from group 1 (ground truth)
#' and \eqn{\tilde{S}(t_k)} is the KM estimate from group 2 (emulated data)
#' at time \eqn{t_k}. Survival probabilities are obtained by evaluating the
#' KM step function at each evaluation point.
#'
#' When \code{subgroup_var} is specified, the aggregate across strata is
#' controlled by \code{method_avg}:
#'
#' \strong{Simple average} (\code{method_avg = "simple"}, default):
#' \deqn{
#'   \mathrm{RMSE}_{\mathrm{avg}} :=
#'   \frac{1}{K} \sum_{x=1}^{K} \mathrm{RMSE}_x
#' }
#'
#' \strong{Weighted average} (\code{method_avg = "weighted"}):
#' \deqn{
#'   \mathrm{RMSE}_{\mathrm{wtd}} :=
#'   \sum_{x=1}^{K} \frac{n_x}{n} \mathrm{RMSE}_x
#' }
#'
#' Note: Lang et al. (2025) defines RMSE for individual strata only and
#' does not specify an aggregation method across strata. The default
#' \code{"simple"} average is used for consistency with
#' \code{\link{ksdist}}.
#'
#' \strong{Subgroup specification via \code{subgroup_var}:}
#'
#' \code{subgroup_var} accepts a character vector where each element is
#' either a single column name or an interaction term written with
#' \code{":"} as separator. The function processes each element
#' independently and stacks all results into a single \code{by_stratum}
#' data frame. The \code{subgroup_spec} column in the output identifies
#' which specification produced each row.
#'
#' \itemize{
#'   \item \code{"SEX"} -- per-stratum results for each level of SEX.
#'   \item \code{"SEX:REGION"} -- per-stratum results for each combination
#'     of SEX and REGION (2 x 3 = 6 strata for the default dummy data).
#'   \item \code{c("SEX", "REGION", "SEX:REGION")} -- all three sets of
#'     results combined in one call.
#'   \item \code{"SEX:REGION:EXTRA"} -- three-way (N-way) interaction;
#'     any number of variables joined by \code{":"} is supported.
#' }
#'
#' When \code{include_overall = TRUE}, the overall population result
#' (i.e., \code{subgroup_spec = "Overall"}) is prepended to \code{by_stratum}
#' and \code{aggregate} regardless of what \code{subgroup_var} specifies.
#' When \code{subgroup_var = NULL}, \code{include_overall} is ignored because
#' the overall result is the only result returned.
#'
#' @param data_group1  Input for group 1 (e.g., original IPD, treated as
#'   ground truth). Accepted formats: a \code{data.frame} or tibble already
#'   in memory, a character path to a \code{.csv} file, or a character path
#'   to an \code{.RData} / \code{.rda} file containing a single object.
#'   Must contain columns specified by \code{time_var}, \code{event_var},
#'   and all variables referenced in \code{subgroup_var}.
#'   See \code{\link{load_ipd}} for details.
#' @param data_group2  Input for group 2 (e.g., emulated IPD). Same accepted
#'   formats as \code{data_group1}.
#' @param time_var     Character. Name of the time-to-event column.
#'   Default \code{"TIME"}.
#' @param event_var    Character. Name of the event indicator column
#'   (1 = event, 0 = censored). Default \code{"EVENT"}.
#' @param tau          Numeric scalar. Upper limit of the evaluation grid.
#'   The \code{n_points} time points are spaced evenly from 0 to \code{tau}.
#'   If \code{NULL} (default), the maximum observed time within each stratum
#'   is used, which corresponds to the definition in Lang et al. (2025).
#' @param n_points     Positive integer. Number of evenly spaced time points
#'   used for RMSE calculation. Default \code{100L}.
#' @param subgroup_var Character vector or \code{NULL}. Each element is
#'   either a single column name (e.g., \code{"SEX"}) or an interaction
#'   term using \code{":"} as separator (e.g., \code{"SEX:REGION"}).
#'   N-way interactions (e.g., \code{"SEX:REGION:VAR3"}) are supported
#'   for any N. Multiple elements are processed independently and their
#'   results are combined. If \code{NULL} (default), RMSE is computed for
#'   the overall population only.
#' @param method_avg   Character. Method for aggregating per-stratum RMSE
#'   values. One of \code{"simple"} (unweighted average, default) or
#'   \code{"weighted"} (sample-size-weighted average). Ignored when
#'   \code{subgroup_var = NULL}.
#' @param include_overall Logical. If \code{TRUE}, the overall population
#'   result (\code{subgroup_spec = "Overall"}) is prepended to the output
#'   in addition to the subgroup results. Ignored when
#'   \code{subgroup_var = NULL} (overall is always returned in that case).
#'   Default \code{FALSE}.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\code{by_stratum}}{A \code{data.frame} with columns
#'     \code{subgroup_spec} (the element of \code{subgroup_var} that
#'     produced this row, or \code{"Overall"}),
#'     \code{stratum} (stratum label within that specification),
#'     \code{n_group1}, \code{n_group2}, \code{n_total}, and
#'     \code{rmse}.}
#'   \item{\code{aggregate}}{A \code{data.frame} with columns
#'     \code{subgroup_spec} and \code{rmse_aggregate}, giving one aggregate
#'     value per element of \code{subgroup_var} (or a single row with
#'     \code{subgroup_spec = "Overall"} when \code{subgroup_var = NULL}).
#'     When \code{include_overall = TRUE}, the overall row is prepended.}
#'   \item{\code{method_avg}}{The aggregation method used.}
#'   \item{\code{tau_input}}{The user-supplied \code{tau} value
#'     (\code{NULL} when per-stratum maximum observation time was used).}
#'   \item{\code{n_points}}{The number of evenly spaced time points used.}
#' }
#'
#' @references
#' Lang, et al. (2025). RESOLVE-IPD: Reconstruction of Survival data with
#' Latent Outcomes from published Evidence for Individual Patient Data
#' meta-analysis.
#'
#' Zhang, et al. (2024). SurvdigitizeR: an algorithm for automated survival
#' curve digitization. BMC Medical Research Methodology, 24:147.
#'
#' @seealso \code{\link{load_ipd}}, \code{\link{km_est}},
#'   \code{\link{nauc}}, \code{\link{ksdist}}
#'
#'
#' @examples
#' # -------------------------------------------------------------------
#' # Prepare dummy data with gen_dummy_data()
#' # -------------------------------------------------------------------
#' ipd <- gen_dummy_data(seed = 42)
#'
#' # Filter to OS endpoint and split by group
#' ipd_os  <- ipd[ipd$TYPE == "OS", ]
#' grp1_os <- ipd_os[ipd_os$GROUP == "Drug",    ]
#' grp2_os <- ipd_os[ipd_os$GROUP == "Control", ]
#'
#' # Filter to PFS endpoint and split by group
#' ipd_pfs  <- ipd[ipd$TYPE == "PFS", ]
#' grp1_pfs <- ipd_pfs[ipd_pfs$GROUP == "Drug",    ]
#' grp2_pfs <- ipd_pfs[ipd_pfs$GROUP == "Control", ]
#'
#' # -------------------------------------------------------------------
#' # 1. Overall population only (subgroup_var = NULL)
#' # -------------------------------------------------------------------
#' res_overall <- rmse(
#'   data_group1 = grp1_os,
#'   data_group2 = grp2_os,
#'   time_var    = "TIME",
#'   event_var   = "EVENT"
#' )
#' res_overall$by_stratum
#' res_overall$aggregate
#'
#' # -------------------------------------------------------------------
#' # 2. Overall population, OS with tau specified
#' # -------------------------------------------------------------------
#' res_tau <- rmse(
#'   data_group1 = grp1_os,
#'   data_group2 = grp2_os,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' res_tau$by_stratum
#'
#' # -------------------------------------------------------------------
#' # 3. Single subgroup: REGION, OS (simple average, default)
#' # -------------------------------------------------------------------
#' res_region <- rmse(
#'   data_group1  = grp1_os,
#'   data_group2  = grp2_os,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 36,
#'   subgroup_var = "REGION"
#' )
#' res_region$by_stratum
#' res_region$aggregate
#'
#' # -------------------------------------------------------------------
#' # 4. REGION with include_overall = TRUE
#' # -------------------------------------------------------------------
#' res_region_ov <- rmse(
#'   data_group1     = grp1_os,
#'   data_group2     = grp2_os,
#'   time_var        = "TIME",
#'   event_var       = "EVENT",
#'   tau             = 36,
#'   subgroup_var    = "REGION",
#'   include_overall = TRUE
#' )
#' res_region_ov$by_stratum
#' res_region_ov$aggregate
#'
#' # -------------------------------------------------------------------
#' # 5. Two-way interaction: SEX x REGION (2 x 3 = 6 strata), OS
#' # -------------------------------------------------------------------
#' res_interaction <- rmse(
#'   data_group1  = grp1_os,
#'   data_group2  = grp2_os,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 36,
#'   subgroup_var = "SEX:REGION"
#' )
#' res_interaction$by_stratum
#' res_interaction$aggregate
#'
#' # -------------------------------------------------------------------
#' # 6. Mixed: SEX, REGION, SEX x REGION, with Overall prepended
#' # -------------------------------------------------------------------
#' res_mixed <- rmse(
#'   data_group1     = grp1_os,
#'   data_group2     = grp2_os,
#'   time_var        = "TIME",
#'   event_var       = "EVENT",
#'   tau             = 36,
#'   subgroup_var    = c("SEX", "REGION", "SEX:REGION"),
#'   include_overall = TRUE
#' )
#' res_mixed$by_stratum
#' res_mixed$aggregate
#'
#' # -------------------------------------------------------------------
#' # 7. Single subgroup: SEX, PFS (weighted average)
#' # -------------------------------------------------------------------
#' res_sex_wtd <- rmse(
#'   data_group1  = grp1_pfs,
#'   data_group2  = grp2_pfs,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 24,
#'   subgroup_var = "SEX",
#'   method_avg   = "weighted"
#' )
#' res_sex_wtd$aggregate
#'
#' # -------------------------------------------------------------------
#' # 8. CSV file input via load_ipd()
#' # -------------------------------------------------------------------
#' tmp1 <- tempfile(fileext = ".csv")
#' tmp2 <- tempfile(fileext = ".csv")
#' write.csv(grp1_os, tmp1, row.names = FALSE)
#' write.csv(grp2_os, tmp2, row.names = FALSE)
#' res_csv <- rmse(
#'   data_group1 = tmp1,
#'   data_group2 = tmp2,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' res_csv$aggregate
#' unlink(c(tmp1, tmp2))
#'
#' # -------------------------------------------------------------------
#' # 9. RData file (multiple objects): use load_ipd() explicitly
#' # -------------------------------------------------------------------
#' tmp_multi <- tempfile(fileext = ".RData")
#' save(grp1_os, grp2_os, file = tmp_multi)
#' d1 <- load_ipd(tmp_multi, rdata_object = "grp1_os")
#' d2 <- load_ipd(tmp_multi, rdata_object = "grp2_os")
#' res_multi <- rmse(
#'   data_group1 = d1,
#'   data_group2 = d2,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' res_multi$aggregate
#' unlink(tmp_multi)
#'
#' @export
rmse <- function(data_group1,
                 data_group2,
                 time_var        = "TIME",
                 event_var       = "EVENT",
                 tau             = NULL,
                 n_points        = 100L,
                 subgroup_var    = NULL,
                 method_avg      = c("simple", "weighted"),
                 include_overall = FALSE) {

  method_avg <- match.arg(method_avg)

  # ---------------------------------------------------------------------------
  # Step 1: Load input data via load_ipd()
  # ---------------------------------------------------------------------------
  data_group1 <- load_ipd(data_group1)
  data_group2 <- load_ipd(data_group2)

  # ---------------------------------------------------------------------------
  # Internal helper: evaluate KM step function at time t
  # km_df is the data.frame(time, surv) returned by km_est()
  # Row i represents S(t) on [time[i], time[i+1])
  # ---------------------------------------------------------------------------
  .km_eval <- function(km_df, t) {
    idx <- findInterval(t, km_df$time)
    if (idx == 0L) return(1.0)
    km_df$surv[idx]
  }

  # ---------------------------------------------------------------------------
  # Internal helper: compute RMSE for one stratum given two KM data.frames
  #
  # Evaluates both KM step functions at n_pts evenly spaced time points
  # from 0 to tau_val, then computes:
  #   RMSE = sqrt( mean( (S1(t_k) - S2(t_k))^2 ) )
  # ---------------------------------------------------------------------------
  .rmse_from_km <- function(km1, km2, tau_val, n_pts) {
    eval_times <- seq(0, tau_val, length.out = n_pts)
    s1 <- sapply(eval_times, function(tt) .km_eval(km1, tt))
    s2 <- sapply(eval_times, function(tt) .km_eval(km2, tt))
    sqrt(mean((s1 - s2)^2))
  }

  # ---------------------------------------------------------------------------
  # Internal helper: build a stratum character vector for one spec string.
  #
  # A spec without ":" is treated as a single column name.
  # A spec with ":" is split on ":" and interaction() is applied, supporting
  # N-way interactions (N >= 2) for any number of variables.
  # Stratum labels are formatted as "v1 / v2 / ..." for interaction terms.
  # ---------------------------------------------------------------------------
  .make_strata_col <- function(data, spec) {
    vars <- trimws(strsplit(spec, ":", fixed = TRUE)[[1]])
    missing_vars <- vars[!vars %in% names(data)]
    if (length(missing_vars) > 0L)
      stop(sprintf(
        "Variable(s) not found in data: %s",
        paste(missing_vars, collapse = ", ")
      ))
    if (length(vars) == 1L) {
      as.character(data[[vars]])
    } else {
      as.character(interaction(data[, vars, drop = FALSE], sep = " / "))
    }
  }

  # ---------------------------------------------------------------------------
  # Internal helper: compute per-stratum RMSE for one specification.
  # Returns a data.frame with columns:
  #   subgroup_spec, stratum, n_group1, n_group2, n_total, rmse
  # ---------------------------------------------------------------------------
  .rmse_one_spec <- function(d1, d2, spec, tau_user, n_pts) {
    if (spec == "Overall") {
      strata_d1 <- rep("Overall", nrow(d1))
      strata_d2 <- rep("Overall", nrow(d2))
    } else {
      strata_d1 <- .make_strata_col(d1, spec)
      strata_d2 <- .make_strata_col(d2, spec)
    }

    strata <- sort(unique(c(strata_d1, strata_d2)))

    rows <- lapply(strata, function(s) {
      sub1 <- d1[strata_d1 == s, ]
      sub2 <- d2[strata_d2 == s, ]
      n1   <- nrow(sub1)
      n2   <- nrow(sub2)

      if (n1 == 0L || n2 == 0L) {
        warning(sprintf(
          "Spec '%s', stratum '%s': empty subset in one group. RMSE set to NA.",
          spec, s
        ))
        return(data.frame(subgroup_spec    = spec,
                          stratum          = as.character(s),
                          n_group1         = n1,
                          n_group2         = n2,
                          n_total          = n1 + n2,
                          rmse             = NA_real_,
                          stringsAsFactors = FALSE))
      }

      # When tau = NULL, use the maximum observed time within the stratum so
      # that each stratum is evaluated over its own full observation range.
      # This makes the per-stratum RMSE consistent with the value displayed
      # in the corresponding kmplot() panel, which receives only stratum data.
      tau_s <- if (is.null(tau_user)) {
        max(c(sub1[[time_var]], sub2[[time_var]]), na.rm = TRUE)
      } else {
        tau_user
      }

      km_res   <- km_est(sub1, sub2,
                         time_var  = time_var,
                         event_var = event_var,
                         tau       = tau_s)
      rmse_val <- .rmse_from_km(km_res$km1, km_res$km2,
                                tau_val = tau_s,
                                n_pts   = n_pts)

      data.frame(subgroup_spec    = spec,
                 stratum          = as.character(s),
                 n_group1         = n1,
                 n_group2         = n2,
                 n_total          = n1 + n2,
                 rmse             = rmse_val,
                 stringsAsFactors = FALSE)
    })

    do.call(rbind, rows)
  }

  # ---------------------------------------------------------------------------
  # Internal helper: compute aggregate RMSE for one specification.
  # Returns a single-row data.frame with columns: subgroup_spec, rmse_aggregate
  # ---------------------------------------------------------------------------
  .aggregate_one_spec <- function(rows, spec) {
    agg_val <- if (method_avg == "weighted") {
      n_grand <- sum(rows$n_total, na.rm = TRUE)
      sum(rows$n_total * rows$rmse, na.rm = TRUE) / n_grand
    } else {
      mean(rows$rmse, na.rm = TRUE)
    }
    data.frame(subgroup_spec   = spec,
               rmse_aggregate  = agg_val,
               stringsAsFactors = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Step 2: Input validation
  # ---------------------------------------------------------------------------
  for (v in c(time_var, event_var)) {
    if (!v %in% names(data_group1))
      stop(sprintf("'%s' not found in data_group1.", v))
    if (!v %in% names(data_group2))
      stop(sprintf("'%s' not found in data_group2.", v))
  }
  if (!is.numeric(n_points) || length(n_points) != 1L || n_points < 2L)
    stop("'n_points' must be a single integer >= 2.")
  if (!is.logical(include_overall) || length(include_overall) != 1L)
    stop("'include_overall' must be a single logical value (TRUE or FALSE).")

  # ---------------------------------------------------------------------------
  # Step 3: Validate tau if explicitly provided
  # When tau = NULL, the upper limit is determined per stratum inside
  # .rmse_one_spec() from the subsetted data.
  # ---------------------------------------------------------------------------
  if (!is.null(tau)) {
    if (!is.numeric(tau) || length(tau) != 1L || tau <= 0)
      stop("'tau' must be a single positive number.")
  }

  # ---------------------------------------------------------------------------
  # Step 4: Build list of specifications to process
  # When subgroup_var = NULL, only "Overall" is processed and
  # include_overall is ignored.
  # When subgroup_var is specified and include_overall = TRUE, "Overall"
  # is prepended to the spec list so the overall result appears first.
  # ---------------------------------------------------------------------------
  if (is.null(subgroup_var)) {
    specs <- "Overall"
  } else {
    specs <- if (include_overall) c("Overall", subgroup_var) else subgroup_var
  }

  n_pts_i <- as.integer(n_points)

  # ---------------------------------------------------------------------------
  # Step 5: Compute per-stratum RMSE for each specification and combine
  # ---------------------------------------------------------------------------
  all_results <- lapply(specs, function(sp) {
    .rmse_one_spec(data_group1, data_group2, sp, tau, n_pts_i)
  })

  by_stratum           <- do.call(rbind, all_results)
  rownames(by_stratum) <- NULL

  # ---------------------------------------------------------------------------
  # Step 6: Compute one aggregate value per specification
  # ---------------------------------------------------------------------------
  agg_rows <- lapply(specs, function(sp) {
    rows <- by_stratum[by_stratum$subgroup_spec == sp, ]
    .aggregate_one_spec(rows, sp)
  })

  aggregate            <- do.call(rbind, agg_rows)
  rownames(aggregate)  <- NULL

  list(
    by_stratum = by_stratum,
    aggregate  = aggregate,
    method_avg = method_avg,
    tau_input  = tau,
    n_points   = n_pts_i
  )
}
