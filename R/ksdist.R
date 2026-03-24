#' Kolmogorov-Smirnov Distance Between Two KM Estimates
#'
#' @description
#' Computes the Kolmogorov-Smirnov (K-S) distance as defined in Definition 2
#' of Zhao et al. (2025, SynthIPD). The K-S distance measures the single
#' largest vertical gap between two Kaplan-Meier survival curves over
#' \eqn{[0, \tau]}. KM estimation is performed internally via
#' \code{\link{km_est}}. Input data can be supplied as a \code{data.frame},
#' tibble, path to a \code{.csv} file, or path to an \code{.RData} /
#' \code{.rda} file containing a single object, via \code{\link{load_ipd}}.
#' For \code{.RData} files containing multiple objects, call
#' \code{\link{load_ipd}} explicitly with the \code{rdata_object} argument
#' before passing the result here.
#'
#' @details
#' The K-S distance for a single stratum is estimated using the discrete
#' estimator. Let the distinct, ordered event times up to \eqn{\tau} across
#' both groups be \eqn{0 = t_{(0)} \le \cdots \le t_{(q)} \le \tau}. Then:
#'
#' \deqn{
#'   \widehat{KS}(\tilde{d}_n \mid d_n) :=
#'   \max_{1 \le j \le q}
#'   \left| \hat{S}(t_{(j)}) - \tilde{S}(t_{(j)}) \right|
#' }
#'
#' When \code{subgroup_var} is specified, the aggregate across strata is
#' controlled by \code{method_avg}:
#'
#' \strong{Simple average} (\code{method_avg = "simple"}, default,
#' per Zhao et al. 2025):
#' \deqn{
#'   \widehat{KS}_{\mathrm{avg}} :=
#'   \frac{1}{K} \sum_{x=1}^{K} \widehat{KS}_x
#' }
#'
#' \strong{Weighted average} (\code{method_avg = "weighted"}):
#' \deqn{
#'   \widehat{KS}_{\mathrm{wtd}} :=
#'   \sum_{x=1}^{K} \frac{n_x}{n} \widehat{KS}_x
#' }
#'
#' Note that Zhao et al. (2025) defines the K-S distance aggregate as a
#' simple average because the K-S distance measures the worst-case scenario
#' and each stratum is considered equally important.
#'
#' The time point at which the maximum gap is achieved is also returned for
#' each stratum, which is useful for annotating KM plots.
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
#' @param data_group1  Input for group 1 (e.g., original IPD). Accepted
#'   formats: a \code{data.frame} or tibble already in memory, a character
#'   path to a \code{.csv} file, or a character path to an \code{.RData} /
#'   \code{.rda} file containing a single object. Must contain columns
#'   specified by \code{time_var}, \code{event_var}, and all variables
#'   referenced in \code{subgroup_var}. See \code{\link{load_ipd}} for
#'   details.
#' @param data_group2  Input for group 2 (e.g., emulated IPD). Same accepted
#'   formats as \code{data_group1}.
#' @param time_var     Character. Name of the time-to-event column.
#'   Default \code{"TIME"}.
#' @param event_var    Character. Name of the event indicator column
#'   (1 = event, 0 = censored). Default \code{"EVENT"}.
#' @param tau          Numeric scalar. Truncation time \eqn{\tau}. The K-S
#'   distance is computed over \eqn{[0, \tau]}. If \code{NULL} (default),
#'   the maximum observed time across both groups is used.
#' @param subgroup_var Character vector or \code{NULL}. Each element is
#'   either a single column name (e.g., \code{"SEX"}) or an interaction
#'   term using \code{":"} as separator (e.g., \code{"SEX:REGION"}).
#'   N-way interactions (e.g., \code{"SEX:REGION:VAR3"}) are supported
#'   for any N. Multiple elements are processed independently and their
#'   results are combined. If \code{NULL} (default), the K-S distance is
#'   computed for the overall population only.
#' @param method_avg   Character. Method for aggregating per-stratum K-S
#'   distance values. One of \code{"simple"} (unweighted average, default,
#'   per Zhao et al. 2025) or \code{"weighted"} (sample-size-weighted
#'   average). Ignored when \code{subgroup_var = NULL}.
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
#'     \code{n_group1}, \code{n_group2}, \code{n_total},
#'     \code{ks} (K-S distance), and
#'     \code{t_max_gap} (time point achieving the maximum gap).}
#'   \item{\code{aggregate}}{A \code{data.frame} with columns
#'     \code{subgroup_spec} and \code{ks_aggregate}, giving one aggregate
#'     value per element of \code{subgroup_var} (or a single row with
#'     \code{subgroup_spec = "Overall"} when \code{subgroup_var = NULL}).
#'     When \code{include_overall = TRUE}, the overall row is prepended.}
#'   \item{\code{method_avg}}{The aggregation method used.}
#'   \item{\code{tau}}{The truncation time actually used.}
#' }
#'
#' @references
#' Zhao, Z. et al. (2025). SynthIPD: Generating synthetic individual
#' patient data from published Kaplan-Meier curves and summary statistics.
#'
#' @seealso \code{\link{load_ipd}}, \code{\link{km_est}},
#'   \code{\link{nauc}}, \code{\link{rmse}}
#'
#' @import dplyr
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
#' res_overall <- ksdist(
#'   data_group1 = grp1_os,
#'   data_group2 = grp2_os,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' res_overall$by_stratum
#' res_overall$aggregate
#'
#' # -------------------------------------------------------------------
#' # 2. Single subgroup: REGION, OS (simple average, default)
#' # -------------------------------------------------------------------
#' res_region <- ksdist(
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
#' # 3. REGION with include_overall = TRUE
#' # -------------------------------------------------------------------
#' res_region_ov <- ksdist(
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
#' # 4. Two-way interaction: SEX x REGION (2 x 3 = 6 strata), OS
#' # -------------------------------------------------------------------
#' res_interaction <- ksdist(
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
#' # 5. Mixed: SEX, REGION, SEX x REGION, with Overall prepended
#' # -------------------------------------------------------------------
#' res_mixed <- ksdist(
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
#' # 6. Single subgroup: SEX, PFS (weighted average)
#' # -------------------------------------------------------------------
#' res_sex_pfs <- ksdist(
#'   data_group1  = grp1_pfs,
#'   data_group2  = grp2_pfs,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 24,
#'   subgroup_var = "SEX",
#'   method_avg   = "weighted"
#' )
#' res_sex_pfs$by_stratum
#'
#' # -------------------------------------------------------------------
#' # 7. CSV file input via load_ipd()
#' # -------------------------------------------------------------------
#' tmp1 <- tempfile(fileext = ".csv")
#' tmp2 <- tempfile(fileext = ".csv")
#' write.csv(grp1_os, tmp1, row.names = FALSE)
#' write.csv(grp2_os, tmp2, row.names = FALSE)
#' res_csv <- ksdist(
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
#' # 8. RData file (multiple objects): use load_ipd() explicitly
#' # -------------------------------------------------------------------
#' tmp_multi <- tempfile(fileext = ".RData")
#' save(grp1_os, grp2_os, file = tmp_multi)
#' d1 <- load_ipd(tmp_multi, rdata_object = "grp1_os")
#' d2 <- load_ipd(tmp_multi, rdata_object = "grp2_os")
#' res_multi <- ksdist(
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
ksdist <- function(data_group1,
                   data_group2,
                   time_var        = "TIME",
                   event_var       = "EVENT",
                   tau             = NULL,
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
  # Internal helper: compute K-S distance for one stratum given two KM
  # data.frames returned by km_est()
  #
  # Uses the discrete estimator from Zhao et al. (2025):
  #   KS = max_{j=1,...,q} |S1(t_j) - S2(t_j)|
  # where t_j are the union of distinct time points within (0, tau_val].
  # t = 0 is excluded because S1(0) = S2(0) = 1 always gives a gap of 0.
  # Also returns t_max_gap: the time point achieving the maximum gap.
  # ---------------------------------------------------------------------------
  .ks_from_km <- function(km1, km2, tau_val) {
    union_times <- sort(unique(c(km1$time, km2$time)))
    union_times <- union_times[union_times > 0 & union_times <= tau_val]

    if (length(union_times) == 0L)
      return(list(ks = 0, t_max_gap = NA_real_))

    abs_diff <- abs(
      sapply(union_times, function(tt) .km_eval(km1, tt)) -
        sapply(union_times, function(tt) .km_eval(km2, tt))
    )
    max_idx <- which.max(abs_diff)
    list(ks = abs_diff[max_idx], t_max_gap = union_times[max_idx])
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
  # Internal helper: compute per-stratum K-S distance for one specification.
  # Returns a data.frame with columns:
  #   subgroup_spec, stratum, n_group1, n_group2, n_total, ks, t_max_gap
  # ---------------------------------------------------------------------------
  .ks_one_spec <- function(d1, d2, spec, tau_val) {
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
          "Spec '%s', stratum '%s': empty subset in one group. KS set to NA.",
          spec, s
        ))
        return(data.frame(subgroup_spec    = spec,
                          stratum          = as.character(s),
                          n_group1         = n1,
                          n_group2         = n2,
                          n_total          = n1 + n2,
                          ks               = NA_real_,
                          t_max_gap        = NA_real_,
                          stringsAsFactors = FALSE))
      }

      km_res <- km_est(sub1, sub2,
                       time_var  = time_var,
                       event_var = event_var,
                       tau       = tau_val)
      res    <- .ks_from_km(km_res$km1, km_res$km2, tau_val = tau_val)

      data.frame(subgroup_spec    = spec,
                 stratum          = as.character(s),
                 n_group1         = n1,
                 n_group2         = n2,
                 n_total          = n1 + n2,
                 ks               = res$ks,
                 t_max_gap        = res$t_max_gap,
                 stringsAsFactors = FALSE)
    })

    do.call(rbind, rows)
  }

  # ---------------------------------------------------------------------------
  # Internal helper: compute aggregate K-S distance for one specification.
  # Returns a single-row data.frame with columns: subgroup_spec, ks_aggregate
  # ---------------------------------------------------------------------------
  .aggregate_one_spec <- function(rows, spec) {
    agg_val <- if (method_avg == "weighted") {
      n_grand <- sum(rows$n_total, na.rm = TRUE)
      sum(rows$n_total * rows$ks, na.rm = TRUE) / n_grand
    } else {
      mean(rows$ks, na.rm = TRUE)
    }
    data.frame(subgroup_spec    = spec,
               ks_aggregate     = agg_val,
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
  if (!is.logical(include_overall) || length(include_overall) != 1L)
    stop("'include_overall' must be a single logical value (TRUE or FALSE).")

  # ---------------------------------------------------------------------------
  # Step 3: Determine tau from the overall data (before subsetting)
  # ---------------------------------------------------------------------------
  t_all <- c(data_group1[[time_var]], data_group2[[time_var]])
  if (is.null(tau)) {
    tau <- max(t_all, na.rm = TRUE)
  } else {
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

  # ---------------------------------------------------------------------------
  # Step 5: Compute per-stratum K-S distance for each specification and combine
  # ---------------------------------------------------------------------------
  all_results <- lapply(specs, function(sp) {
    .ks_one_spec(data_group1, data_group2, sp, tau)
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
    tau        = tau
  )
}
