#' Compute NAUC, K-S Distance, RMSE, and Descriptive Statistics by Subgroup
#'
#' @description
#' Computes three evaluation metrics -- NAUC, Kolmogorov-Smirnov (K-S)
#' distance, and RMSE -- together with descriptive survival statistics
#' (sample size, event count, median survival time with 95\% CI) for the
#' overall population and for each level of one or more subgroup variables.
#' Results are returned as an object of class \code{EvalMetrics}, which
#' inherits from \code{data.frame} and has a dedicated \code{print()} method
#' for formatted console output.
#'
#' @details
#' Each metric is computed by calling \code{\link{nauc}},
#' \code{\link{ksdist}}, and \code{\link{rmse}} respectively with
#' \code{include_overall = TRUE}, so the overall population result is always
#' prepended to each subgroup block.
#'
#' NAUC and K-S distance are computed over \eqn{[0, \tau]}. RMSE is always
#' computed over the full observation range of each stratum (\code{tau = NULL}
#' internally), following Lang et al. (2025).
#'
#' Descriptive statistics (sample size, events, median survival time and its
#' CI) are obtained from \code{\link{km_est}}, which implements the Greenwood
#' variance formula with a log transformation for the CI -- the same approach
#' as \code{survival::survfit} with \code{conf.type = "log"}.
#'
#' \strong{Subgroup specification via \code{subgroup_vars}:}
#'
#' Each element of \code{subgroup_vars} follows the same \code{":"} separator
#' convention as \code{\link{nauc}}, \code{\link{ksdist}}, and
#' \code{\link{rmse}}: a single column name (e.g., \code{"SEX"}) or an
#' interaction term (e.g., \code{"SEX:REGION"}). N-way interactions are
#' supported for any N.
#'
#' \strong{S3 class \code{EvalMetrics}:}
#'
#' The returned object has class \code{c("EvalMetrics", "data.frame")} and
#' can be used as a plain \code{data.frame} in all downstream operations.
#' Calling \code{print()} on it renders a formatted table grouped by
#' \code{subgroup_spec}, with descriptive statistics and metric values
#' side by side.
#'
#' @param data_group1      Input for group 1 (e.g., original IPD). Accepted
#'   formats: a \code{data.frame} or tibble already in memory, a character
#'   path to a \code{.csv} file, or a character path to an \code{.RData} /
#'   \code{.rda} file containing a single object. Must contain columns
#'   specified by \code{time_var}, \code{event_var}, and all variables
#'   referenced in \code{subgroup_vars}. See \code{\link{load_ipd}} for
#'   details.
#' @param data_group2      Input for group 2 (e.g., emulated IPD). Same
#'   accepted formats as \code{data_group1}.
#' @param time_var         Character. Name of the time-to-event column.
#'   Default \code{"TIME"}.
#' @param event_var        Character. Name of the event indicator column
#'   (1 = event, 0 = censored). Default \code{"EVENT"}.
#' @param tau              Numeric scalar. Truncation time \eqn{\tau} used
#'   for NAUC and K-S distance. If \code{NULL} (default), the maximum
#'   observed time across both groups is used. RMSE is always computed over
#'   the full observation range regardless of this value.
#' @param subgroup_vars    Character vector or \code{NULL}. Each element is
#'   either a single column name (e.g., \code{"SEX"}) or an interaction term
#'   using \code{":"} as separator (e.g., \code{"SEX:REGION"}). N-way
#'   interactions are supported for any N. If \code{NULL} (default), only the
#'   overall population is computed.
#' @param n_points         Positive integer. Number of evenly spaced time
#'   points used for RMSE calculation. Default \code{100L}.
#' @param conf_level       Numeric in \code{(0, 1)}. Confidence level for
#'   the median survival time CI. Default \code{0.95}. Passed to
#'   \code{\link{km_est}}.
#' @param label_group1     Character. Display label for group 1 used in the
#'   \code{print()} output. Default \code{"g1"}.
#' @param label_group2     Character. Display label for group 2 used in the
#'   \code{print()} output. Default \code{"g2"}.
#' @param nauc_method_avg  Character. Aggregation method for NAUC across
#'   strata. One of \code{"weighted"} (default, per Zhao et al. 2025) or
#'   \code{"simple"}. Passed to \code{\link{nauc}}.
#' @param ks_method_avg    Character. Aggregation method for K-S distance
#'   across strata. One of \code{"simple"} (default, per Zhao et al. 2025)
#'   or \code{"weighted"}. Passed to \code{\link{ksdist}}.
#' @param rmse_method_avg  Character. Aggregation method for RMSE across
#'   strata. One of \code{"simple"} (default) or \code{"weighted"}.
#'   Passed to \code{\link{rmse}}.
#'
#' @return An object of class \code{c("EvalMetrics", "data.frame")} with one
#'   row per stratum per subgroup specification. Columns:
#' \describe{
#'   \item{\code{subgroup_spec}}{Subgroup specification string, or
#'     \code{"Overall"} for the overall population block.}
#'   \item{\code{stratum}}{Level of the subgroup variable, or
#'     \code{"Overall"} for the overall population row.}
#'   \item{\code{n_group1}}{Sample size of group 1 in that stratum.}
#'   \item{\code{n_group2}}{Sample size of group 2 in that stratum.}
#'   \item{\code{events_group1}}{Number of events in group 1.}
#'   \item{\code{events_group2}}{Number of events in group 2.}
#'   \item{\code{median_group1}}{Median survival time for group 1.}
#'   \item{\code{median_ci_group1}}{95\% CI for median of group 1 as a
#'     character string \code{"[lower, upper]"}.}
#'   \item{\code{median_group2}}{Median survival time for group 2.}
#'   \item{\code{median_ci_group2}}{95\% CI for median of group 2 as a
#'     character string \code{"[lower, upper]"}.}
#'   \item{\code{nauc}}{NAUC value for that stratum.}
#'   \item{\code{ks}}{K-S distance value for that stratum.}
#'   \item{\code{t_max_gap}}{Time point at which the maximum K-S gap is
#'     achieved.}
#'   \item{\code{rmse}}{RMSE value for that stratum.}
#'   \item{\code{nauc_aggregate}}{Aggregate NAUC across strata within the
#'     subgroup specification. \code{NA} for the \code{"Overall"} block.}
#'   \item{\code{ks_aggregate}}{Aggregate K-S distance. \code{NA} for the
#'     \code{"Overall"} block.}
#'   \item{\code{rmse_aggregate}}{Aggregate RMSE. \code{NA} for the
#'     \code{"Overall"} block.}
#' }
#'
#' @references
#' Zhao, Z. et al. (2025). SynthIPD: Generating synthetic individual
#' patient data from published Kaplan-Meier curves and summary statistics.
#'
#' Lang, et al. (2025). RESOLVE-IPD: Reconstruction of Survival data with
#' Latent Outcomes from published Evidence for Individual Patient Data
#' meta-analysis.
#'
#' @seealso \code{\link{load_ipd}}, \code{\link{km_est}}, \code{\link{nauc}},
#'   \code{\link{ksdist}}, \code{\link{rmse}}
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
#' # 1. Overall only
#' # -------------------------------------------------------------------
#' res <- summarytable(
#'   data_group1  = grp1_os,
#'   data_group2  = grp2_os,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 36,
#'   label_group1 = "Original",
#'   label_group2 = "Emulated"
#' )
#' print(res)
#'
#' # -------------------------------------------------------------------
#' # 2. Overall + single subgroup: REGION, OS
#' # -------------------------------------------------------------------
#' res_region <- summarytable(
#'   data_group1   = grp1_os,
#'   data_group2   = grp2_os,
#'   time_var      = "TIME",
#'   event_var     = "EVENT",
#'   tau           = 36,
#'   subgroup_vars = "REGION",
#'   label_group1  = "Original",
#'   label_group2  = "Emulated"
#' )
#' print(res_region)
#'
#' # -------------------------------------------------------------------
#' # 3. Overall + multiple subgroup variables, OS
#' # -------------------------------------------------------------------
#' res_os <- summarytable(
#'   data_group1   = grp1_os,
#'   data_group2   = grp2_os,
#'   time_var      = "TIME",
#'   event_var     = "EVENT",
#'   tau           = 36,
#'   subgroup_vars = c("SEX", "REGION"),
#'   label_group1  = "Original",
#'   label_group2  = "Emulated"
#' )
#' print(res_os)
#'
#' # -------------------------------------------------------------------
#' # 4. Two-way interaction: SEX x REGION, OS
#' # -------------------------------------------------------------------
#' res_inter <- summarytable(
#'   data_group1   = grp1_os,
#'   data_group2   = grp2_os,
#'   time_var      = "TIME",
#'   event_var     = "EVENT",
#'   tau           = 36,
#'   subgroup_vars = "SEX:REGION",
#'   label_group1  = "Original",
#'   label_group2  = "Emulated"
#' )
#' print(res_inter)
#'
#' # -------------------------------------------------------------------
#' # 5. Mixed specification, PFS, weighted RMSE
#' # -------------------------------------------------------------------
#' res_pfs <- summarytable(
#'   data_group1      = grp1_pfs,
#'   data_group2      = grp2_pfs,
#'   time_var         = "TIME",
#'   event_var        = "EVENT",
#'   tau              = 24,
#'   subgroup_vars    = c("SEX", "REGION", "SEX:REGION"),
#'   label_group1     = "Original",
#'   label_group2     = "Emulated",
#'   rmse_method_avg  = "weighted"
#' )
#' print(res_pfs)
#'
#' # Access as a plain data.frame
#' as.data.frame(res_pfs)
#'
#' @export
summarytable <- function(data_group1,
                         data_group2,
                         time_var        = "TIME",
                         event_var       = "EVENT",
                         tau             = NULL,
                         subgroup_vars   = NULL,
                         n_points        = 100L,
                         conf_level      = 0.95,
                         label_group1    = "g1",
                         label_group2    = "g2",
                         nauc_method_avg = c("weighted", "simple"),
                         ks_method_avg   = c("simple", "weighted"),
                         rmse_method_avg = c("simple", "weighted")) {

  nauc_method_avg <- match.arg(nauc_method_avg)
  ks_method_avg   <- match.arg(ks_method_avg)
  rmse_method_avg <- match.arg(rmse_method_avg)

  # ---------------------------------------------------------------------------
  # Step 1: Load input data via load_ipd()
  # ---------------------------------------------------------------------------
  data_group1 <- load_ipd(data_group1)
  data_group2 <- load_ipd(data_group2)

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
  if (!is.null(tau)) {
    if (!is.numeric(tau) || length(tau) != 1L || tau <= 0)
      stop("'tau' must be a single positive number.")
  }

  # ---------------------------------------------------------------------------
  # Step 3: Determine tau for NAUC and K-S distance
  # RMSE always uses tau = NULL (full range per Lang et al. 2025)
  # ---------------------------------------------------------------------------
  t_all       <- c(data_group1[[time_var]], data_group2[[time_var]])
  tau_nauc_ks <- if (is.null(tau)) max(t_all, na.rm = TRUE) else tau

  # ---------------------------------------------------------------------------
  # Internal helper: build a stratum character vector for one spec string.
  # ---------------------------------------------------------------------------
  .make_strata_col <- function(data, spec) {
    vars         <- trimws(strsplit(spec, ":", fixed = TRUE)[[1]])
    missing_vars <- vars[!vars %in% names(data)]
    if (length(missing_vars) > 0L)
      stop(sprintf("Variable(s) not found in data: %s",
                   paste(missing_vars, collapse = ", ")))
    if (length(vars) == 1L) {
      as.character(data[[vars]])
    } else {
      as.character(interaction(data[, vars, drop = FALSE], sep = " / "))
    }
  }

  # ---------------------------------------------------------------------------
  # Internal helper: extract descriptive stats from km_est() for one stratum.
  # Returns a one-row data.frame with desc columns for both groups.
  # ---------------------------------------------------------------------------
  .desc_row <- function(d1, d2, subgroup_spec, stratum) {
    km <- km_est(d1, d2,
                 time_var   = time_var,
                 event_var  = event_var,
                 conf_level = conf_level)

    # Format CI as "[lower, upper]" with NA handled gracefully
    .fmt_ci <- function(lo, hi, dg = 2L) {
      lo_s <- if (is.na(lo)) "NA" else formatC(lo, format = "f", digits = dg)
      hi_s <- if (is.na(hi)) "NA" else formatC(hi, format = "f", digits = dg)
      sprintf("[%s, %s]", lo_s, hi_s)
    }

    data.frame(
      subgroup_spec    = subgroup_spec,
      stratum          = stratum,
      n_group1         = km$desc1$n,
      n_group2         = km$desc2$n,
      events_group1    = km$desc1$n_events,
      events_group2    = km$desc2$n_events,
      median_group1    = km$desc1$median,
      median_ci_group1 = .fmt_ci(km$desc1$ci_lower, km$desc1$ci_upper),
      median_group2    = km$desc2$median,
      median_ci_group2 = .fmt_ci(km$desc2$ci_lower, km$desc2$ci_upper),
      stringsAsFactors = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Internal helper: collect desc rows for all strata of one specification.
  # ---------------------------------------------------------------------------
  .desc_one_spec <- function(spec) {
    if (spec == "Overall") {
      return(.desc_row(data_group1, data_group2, "Overall", "Overall"))
    }

    strata_d1 <- .make_strata_col(data_group1, spec)
    strata_d2 <- .make_strata_col(data_group2, spec)
    strata    <- sort(unique(c(strata_d1, strata_d2)))

    rows <- lapply(strata, function(s) {
      sub1 <- data_group1[strata_d1 == s, ]
      sub2 <- data_group2[strata_d2 == s, ]
      if (nrow(sub1) == 0L || nrow(sub2) == 0L) {
        warning(sprintf(
          "Spec '%s', stratum '%s': empty subset in one group. Skipped.", s))
        return(NULL)
      }
      .desc_row(sub1, sub2, spec, s)
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }

  # ---------------------------------------------------------------------------
  # Step 4: Compute metrics and descriptive statistics for each specification
  # ---------------------------------------------------------------------------
  specs <- if (is.null(subgroup_vars)) character(0) else subgroup_vars

  # -- Overall block: metrics --
  res_nauc_ov <- nauc(
    data_group1 = data_group1, data_group2 = data_group2,
    time_var    = time_var,    event_var   = event_var,
    tau         = tau_nauc_ks, method_avg  = nauc_method_avg
  )
  res_ks_ov <- ksdist(
    data_group1 = data_group1, data_group2 = data_group2,
    time_var    = time_var,    event_var   = event_var,
    tau         = tau_nauc_ks, method_avg  = ks_method_avg
  )
  res_rmse_ov <- rmse(
    data_group1 = data_group1, data_group2 = data_group2,
    time_var    = time_var,    event_var   = event_var,
    tau         = NULL,        n_points    = as.integer(n_points),
    method_avg  = rmse_method_avg
  )

  # -- Overall block: descriptive stats --
  desc_ov <- .desc_one_spec("Overall")

  overall_block <- .merge_results(
    res_nauc   = res_nauc_ov,
    res_ks     = res_ks_ov,
    res_rmse   = res_rmse_ov,
    desc_df    = desc_ov,
    is_overall = TRUE
  )

  # -- Per-subgroup blocks --
  subgroup_blocks <- lapply(specs, function(sp) {
    res_nauc_sg <- nauc(
      data_group1  = data_group1, data_group2  = data_group2,
      time_var     = time_var,    event_var    = event_var,
      tau          = tau_nauc_ks, subgroup_var = sp,
      method_avg   = nauc_method_avg
    )
    res_ks_sg <- ksdist(
      data_group1  = data_group1, data_group2  = data_group2,
      time_var     = time_var,    event_var    = event_var,
      tau          = tau_nauc_ks, subgroup_var = sp,
      method_avg   = ks_method_avg
    )
    res_rmse_sg <- rmse(
      data_group1  = data_group1, data_group2  = data_group2,
      time_var     = time_var,    event_var    = event_var,
      tau          = NULL,        n_points     = as.integer(n_points),
      subgroup_var = sp,          method_avg   = rmse_method_avg
    )
    desc_sg <- .desc_one_spec(sp)

    .merge_results(
      res_nauc   = res_nauc_sg,
      res_ks     = res_ks_sg,
      res_rmse   = res_rmse_sg,
      desc_df    = desc_sg,
      is_overall = FALSE
    )
  })

  # ---------------------------------------------------------------------------
  # Step 5: Combine all blocks and reorder columns
  # ---------------------------------------------------------------------------
  out <- do.call(rbind, c(list(overall_block), subgroup_blocks))
  rownames(out) <- NULL

  col_order <- c("subgroup_spec", "stratum",
                 "n_group1", "n_group2",
                 "events_group1", "events_group2",
                 "median_group1", "median_ci_group1",
                 "median_group2", "median_ci_group2",
                 "nauc", "ks", "t_max_gap", "rmse",
                 "nauc_aggregate", "ks_aggregate", "rmse_aggregate")
  out <- out[, col_order]

  # ---------------------------------------------------------------------------
  # Step 6: Attach metadata and S3 class
  # ---------------------------------------------------------------------------
  attr(out, "tau")             <- tau_nauc_ks
  attr(out, "nauc_method_avg") <- nauc_method_avg
  attr(out, "ks_method_avg")   <- ks_method_avg
  attr(out, "rmse_method_avg") <- rmse_method_avg
  attr(out, "n_points")        <- as.integer(n_points)
  attr(out, "conf_level")      <- conf_level
  attr(out, "time_var")        <- time_var
  attr(out, "event_var")       <- event_var
  attr(out, "label_group1")    <- label_group1
  attr(out, "label_group2")    <- label_group2

  class(out) <- c("EvalMetrics", "data.frame")
  out
}

# =============================================================================
# Internal helper: merge by_stratum outputs from nauc(), ksdist(), rmse()
# and descriptive statistics into a single data.frame with aggregate columns.
# Not exported.
# =============================================================================
.merge_results <- function(res_nauc, res_ks, res_rmse, desc_df,
                           is_overall) {

  # Metric columns from by_stratum
  df_nauc <- res_nauc$by_stratum[, c("subgroup_spec", "stratum", "nauc")]
  df_ks   <- res_ks$by_stratum[,   c("subgroup_spec", "stratum",
                                     "ks", "t_max_gap")]
  df_rmse <- res_rmse$by_stratum[, c("subgroup_spec", "stratum", "rmse")]

  # Merge metrics on subgroup_spec + stratum
  df <- merge(df_nauc, df_ks,   by = c("subgroup_spec", "stratum"),
              sort = FALSE)
  df <- merge(df,      df_rmse, by = c("subgroup_spec", "stratum"),
              sort = FALSE)

  # Merge descriptive stats
  df <- merge(desc_df, df, by = c("subgroup_spec", "stratum"), sort = FALSE)

  # Aggregate columns: NA for Overall (single stratum, no aggregation needed)
  if (is_overall) {
    df$nauc_aggregate <- NA_real_
    df$ks_aggregate   <- NA_real_
    df$rmse_aggregate <- NA_real_
  } else {
    agg_nauc <- res_nauc$aggregate[, c("subgroup_spec", "nauc_aggregate")]
    agg_ks   <- res_ks$aggregate[,   c("subgroup_spec", "ks_aggregate")]
    agg_rmse <- res_rmse$aggregate[, c("subgroup_spec", "rmse_aggregate")]

    df <- merge(df, agg_nauc, by = "subgroup_spec", sort = FALSE)
    df <- merge(df, agg_ks,   by = "subgroup_spec", sort = FALSE)
    df <- merge(df, agg_rmse, by = "subgroup_spec", sort = FALSE)
  }

  df
}

# =============================================================================
#' Print Method for \code{EvalMetrics} Objects
#'
#' @description
#' Prints a formatted summary of an \code{EvalMetrics} object returned by
#' \code{\link{summarytable}}. Results are grouped by \code{subgroup_spec}
#' and displayed in two sections per block: descriptive statistics (N, events,
#' median survival time with CI) followed by evaluation metrics (NAUC, KS,
#' RMSE). The line width is computed dynamically from the content so that
#' separators always fit the output exactly.
#'
#' When \code{kable = TRUE}, the table is formatted using
#' \code{knitr::kable()}. For \code{kable_format = "html"}, additional
#' styling and row grouping are applied via \code{kableExtra} if available.
#'
#' @param x            An object of class \code{EvalMetrics}.
#' @param digits       Non-negative integer. Number of decimal places for
#'   NAUC, K-S distance, and RMSE values. Default \code{4}.
#' @param digits_surv  Non-negative integer. Number of decimal places for
#'   median survival time and CI values. Default \code{2}.
#' @param kable        Logical. If \code{TRUE}, format output using
#'   \code{knitr::kable()}. If \code{FALSE} (default), display as plain
#'   text in the console.
#' @param kable_format Character. Format for \code{kable} output when
#'   \code{kable = TRUE}. One of \code{"pipe"} (default, Markdown),
#'   \code{"simple"}, \code{"latex"}, or \code{"html"}. When
#'   \code{"html"} and \code{kableExtra} is available, additional styling
#'   and subgroup row-grouping are applied via \code{kableExtra::pack_rows()}.
#' @param ...          Currently ignored.
#'
#' @return \code{x} invisibly.
#'
#' @importFrom knitr kable
#'
#' @export
print.EvalMetrics <- function(x, digits = 4L, digits_surv = 2L,
                              kable = FALSE, kable_format = "pipe", ...) {

  kable_format <- match.arg(kable_format,
                            c("pipe", "simple", "latex", "html"))

  digits      <- as.integer(digits)
  digits_surv <- as.integer(digits_surv)
  tau_val     <- attr(x, "tau")
  time_var    <- attr(x, "time_var")
  ev_var      <- attr(x, "event_var")
  n_meth      <- attr(x, "nauc_method_avg")
  k_meth      <- attr(x, "ks_method_avg")
  r_meth      <- attr(x, "rmse_method_avg")
  n_pts       <- attr(x, "n_points")
  conf_lv     <- attr(x, "conf_level")
  lbl1        <- attr(x, "label_group1")
  lbl2        <- attr(x, "label_group2")

  # ---------------------------------------------------------------------------
  # Build header lines
  # ---------------------------------------------------------------------------
  tau_str  <- if (is.null(tau_val)) "max observed" else
    formatC(tau_val, format = "f", digits = 1L)
  agg_str  <- sprintf("NAUC = %s | KS = %s | RMSE = %s",
                      n_meth, k_meth, r_meth)
  rmse_str <- sprintf("full observation range (#points = %d)", n_pts)
  ci_str   <- sprintf("%d%% CI (Greenwood + log transformation)",
                      as.integer(conf_lv * 100))

  hdr_lines <- c(
    sprintf(" Time variable  : %s", time_var),
    sprintf(" Event variable : %s", ev_var),
    sprintf(" tau (NAUC/KS)  : %s", tau_str),
    sprintf(" RMSE range     : %s", rmse_str),
    sprintf(" Median CI      : %s", ci_str),
    sprintf(" Aggregation    : %s", agg_str)
  )

  specs <- unique(x$subgroup_spec)

  # ---------------------------------------------------------------------------
  # Plain text output
  # ---------------------------------------------------------------------------
  if (!kable) {

    # Column widths
    stratum_w    <- max(nchar(x$stratum), nchar("[Aggregate]"),
                        nchar("Stratum")) + 1L
    lbl_w        <- max(nchar(lbl1), nchar(lbl2), nchar("Group")) + 1L
    n_w          <- 6L
    ev_w         <- 6L
    median_w     <- max(digits_surv + 4L, nchar("Median")) + 1L
    ci_w         <- max(2L * (digits_surv + 4L) + 4L, nchar("CI")) + 1L
    metric_w     <- max(digits + 2L, nchar("NAUC"), nchar("KS"),
                        nchar("RMSE")) + 1L

    col_total_w <- 1L + stratum_w + 1L + lbl_w +
      1L + n_w + 1L + ev_w +
      1L + median_w + 1L + ci_w +
      2L + metric_w + 1L + metric_w + 1L + metric_w

    total_w  <- max(col_total_w, max(nchar(hdr_lines)) + 1L)
    sep_eq   <- paste(rep("=", total_w), collapse = "")
    sep_dash <- paste(rep("-", total_w), collapse = "")

    # Format strings
    fmt_hdr <- sprintf(
      " %%-%ds %%-%ds %%%ds %%%ds %%%ds %%-%ds  %%%ds %%%ds %%%ds\n",
      stratum_w, lbl_w, n_w, ev_w, median_w, ci_w,
      metric_w, metric_w, metric_w)

    fmt_row <- sprintf(
      "  %%-%ds %%-%ds %%%dd %%%dd %%%d.%df %%-%ds  %%%d.%df %%%d.%df %%%d.%df\n",
      stratum_w - 1L, lbl_w - 1L, n_w, ev_w,
      median_w, digits_surv, ci_w,
      metric_w, digits, metric_w, digits, metric_w, digits)

    fmt_agg <- sprintf(
      "  %%-%ds %%-%ds %%%ds %%%ds %%%ds %%-%ds  %%%d.%df %%%d.%df %%%d.%df\n",
      stratum_w - 1L, lbl_w - 1L, n_w, ev_w, median_w, ci_w,
      metric_w, digits, metric_w, digits, metric_w, digits)

    cat(sep_eq, "\n")
    cat(" EvalEmulateSurv: Emulation Evaluation Metrics\n")
    cat(sep_eq, "\n")
    for (hl in hdr_lines) cat(hl, "\n")
    cat(sep_dash, "\n")
    cat(sprintf(fmt_hdr,
                "Stratum", "Group", "N", "Events",
                "Median", paste0(as.integer(conf_lv * 100), "% CI"),
                "NAUC", "KS", "RMSE"))
    cat(sep_dash, "\n")

    for (sp in specs) {
      rows <- x[x$subgroup_spec == sp, ]
      cat(sprintf(" [%s]\n", sp))

      for (i in seq_len(nrow(rows))) {
        r       <- rows[i, ]
        med1    <- if (is.na(r$median_group1)) NA_real_ else r$median_group1
        med2    <- if (is.na(r$median_group2)) NA_real_ else r$median_group2

        # Group 1 row
        cat(sprintf(fmt_row,
                    r$stratum, lbl1,
                    r$n_group1, r$events_group1,
                    med1, r$median_ci_group1,
                    r$nauc, r$ks, r$rmse))
        # Group 2 row (stratum and metrics blank)
        cat(sprintf(fmt_row,
                    "", lbl2,
                    r$n_group2, r$events_group2,
                    med2, r$median_ci_group2,
                    r$nauc, r$ks, r$rmse))
      }

      if (sp != "Overall" && !is.na(rows$nauc_aggregate[1])) {
        cat(sprintf(fmt_agg,
                    "[Aggregate]", "", "", "", "", "",
                    rows$nauc_aggregate[1],
                    rows$ks_aggregate[1],
                    rows$rmse_aggregate[1]))
      }
      cat("\n")
    }
    cat(sep_eq, "\n")

    # ---------------------------------------------------------------------------
    # kable output
    # ---------------------------------------------------------------------------
  } else {

    fmt_surv <- function(v) {
      if (is.na(v)) "NA" else formatC(v, format = "f", digits = digits_surv)
    }
    fmt_metric <- function(v) formatC(v, format = "f", digits = digits)

    rows_list <- lapply(specs, function(sp) {
      rows <- x[x$subgroup_spec == sp, ]

      # Two rows per stratum (group1 / group2), metric repeated
      df_str <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
        r <- rows[i, ]
        rbind(
          data.frame(
            Stratum  = r$stratum, Group = lbl1,
            N        = r$n_group1, Events = r$events_group1,
            Median   = fmt_surv(r$median_group1),
            CI       = r$median_ci_group1,
            NAUC     = fmt_metric(r$nauc),
            KS       = fmt_metric(r$ks),
            RMSE     = fmt_metric(r$rmse),
            check.names = FALSE, stringsAsFactors = FALSE),
          data.frame(
            Stratum  = "", Group = lbl2,
            N        = r$n_group2, Events = r$events_group2,
            Median   = fmt_surv(r$median_group2),
            CI       = r$median_ci_group2,
            NAUC     = fmt_metric(r$nauc),
            KS       = fmt_metric(r$ks),
            RMSE     = fmt_metric(r$rmse),
            check.names = FALSE, stringsAsFactors = FALSE)
        )
      }))

      if (sp != "Overall" && !is.na(rows$nauc_aggregate[1])) {
        agg_row <- data.frame(
          Stratum = "[Aggregate]", Group = "",
          N = "", Events = "", Median = "", CI = "",
          NAUC  = fmt_metric(rows$nauc_aggregate[1]),
          KS    = fmt_metric(rows$ks_aggregate[1]),
          RMSE  = fmt_metric(rows$rmse_aggregate[1]),
          check.names = FALSE, stringsAsFactors = FALSE)
        df_str <- rbind(df_str, agg_row)
      }
      df_str
    })

    kable_df           <- do.call(rbind, rows_list)
    rownames(kable_df) <- NULL

    ci_col_name <- sprintf("%d%% CI", as.integer(conf_lv * 100))
    names(kable_df)[names(kable_df) == "CI"] <- ci_col_name

    tbl <- knitr::kable(
      kable_df,
      format  = kable_format,
      escape  = FALSE,
      align   = c("l", "l", "r", "r", "r", "l", "r", "r", "r"),
      caption = sprintf(
        paste0("Emulation Evaluation Metrics | ",
               "time: %s | event: %s | tau: %s | ",
               "RMSE: full range (#pts=%d) | ",
               "Median CI: %d%% Greenwood+log | ",
               "Agg: NAUC=%s, KS=%s, RMSE=%s"),
        time_var, ev_var, tau_str, n_pts,
        as.integer(conf_lv * 100),
        n_meth, k_meth, r_meth
      )
    )

    if (kable_format == "html" &&
        requireNamespace("kableExtra", quietly = TRUE)) {

      tbl <- kableExtra::kable_styling(
        tbl,
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width        = FALSE,
        position          = "center"
      )

      row_offset <- 0L
      for (sp in specs) {
        rows    <- x[x$subgroup_spec == sp, ]
        n_rows  <- nrow(rows)
        has_agg <- sp != "Overall" && !is.na(rows$nauc_aggregate[1])
        # 2 rows per stratum (group1 + group2) + optional aggregate
        n_block <- n_rows * 2L + if (has_agg) 1L else 0L

        tbl <- kableExtra::pack_rows(
          tbl, sp,
          row_offset + 1L,
          row_offset + n_block,
          bold          = TRUE,
          label_row_css = "border-top: 2px solid #555; font-weight: bold;"
        )
        row_offset <- row_offset + n_block
      }
    }

    print(tbl)
  }

  invisible(x)
}
