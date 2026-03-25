#' Compute NAUC, K-S Distance, and RMSE by Subgroup
#'
#' @description
#' Computes three evaluation metrics -- NAUC, Kolmogorov-Smirnov (K-S)
#' distance, and RMSE -- for the overall population and for each level of
#' one or more subgroup variables. Results are returned as an object of
#' class \code{EvalMetrics}, which inherits from \code{data.frame} and has
#' a dedicated \code{print()} method for formatted console output.
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
#' \code{subgroup_spec}, with benchmark indicators for each metric.
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
#'   \item{\code{n_total}}{Combined sample size in that stratum.}
#'   \item{\code{nauc}}{NAUC value for that stratum.}
#'   \item{\code{ks}}{K-S distance value for that stratum.}
#'   \item{\code{t_max_gap}}{Time point at which the maximum K-S gap is
#'     achieved.}
#'   \item{\code{rmse}}{RMSE value for that stratum.}
#'   \item{\code{nauc_aggregate}}{Aggregate NAUC across strata within the
#'     subgroup specification, computed using \code{nauc_method_avg}.
#'     \code{NA} for the \code{"Overall"} block.}
#'   \item{\code{ks_aggregate}}{Aggregate K-S distance across strata,
#'     computed using \code{ks_method_avg}. \code{NA} for the
#'     \code{"Overall"} block.}
#'   \item{\code{rmse_aggregate}}{Aggregate RMSE across strata, computed
#'     using \code{rmse_method_avg}. \code{NA} for the \code{"Overall"}
#'     block.}
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
#' # 1. Overall only
#' # -------------------------------------------------------------------
#' res <- summarytable(
#'   data_group1 = grp1_os,
#'   data_group2 = grp2_os,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
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
#'   subgroup_vars = "REGION"
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
#'   subgroup_vars = c("SEX", "REGION")
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
#'   subgroup_vars = "SEX:REGION"
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
                         nauc_method_avg = c("weighted", "simple"),
                         ks_method_avg   = c("simple", "weighted"),
                         rmse_method_avg = c("simple",
                                             "weighted")) {

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
  # Step 4: Compute all three metrics for each subgroup specification.
  #
  # nauc(), ksdist(), rmse() are each called once per subgroup_vars element
  # with include_overall = TRUE, so the Overall block is always prepended.
  # When subgroup_vars = NULL, only the Overall block is produced.
  # ---------------------------------------------------------------------------
  specs <- if (is.null(subgroup_vars)) character(0) else subgroup_vars

  # -- Overall block --
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

  overall_block <- .merge_results(
    res_nauc  = res_nauc_ov,
    res_ks    = res_ks_ov,
    res_rmse  = res_rmse_ov,
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
    .merge_results(
      res_nauc   = res_nauc_sg,
      res_ks     = res_ks_sg,
      res_rmse   = res_rmse_sg,
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
  attr(out, "time_var")        <- time_var
  attr(out, "event_var")       <- event_var

  class(out) <- c("EvalMetrics", "data.frame")
  out
}

# =============================================================================
# Internal helper: merge by_stratum outputs from nauc(), ksdist(), rmse()
# into a single data.frame with aggregate columns.
# Not exported.
# =============================================================================
.merge_results <- function(res_nauc, res_ks, res_rmse, is_overall) {

  # Extract by_stratum columns needed
  df_nauc <- res_nauc$by_stratum[, c("subgroup_spec", "stratum",
                                     "n_group1", "n_group2",
                                     "nauc")]
  df_ks   <- res_ks$by_stratum[,   c("subgroup_spec", "stratum",
                                     "ks", "t_max_gap")]
  df_rmse <- res_rmse$by_stratum[, c("subgroup_spec", "stratum", "rmse")]

  # Merge on subgroup_spec + stratum
  df <- merge(df_nauc, df_ks,   by = c("subgroup_spec", "stratum"),
              sort = FALSE)
  df <- merge(df,      df_rmse, by = c("subgroup_spec", "stratum"),
              sort = FALSE)

  # Aggregate columns: NA for Overall (single stratum, no aggregation needed)
  if (is_overall) {
    df$nauc_aggregate <- NA_real_
    df$ks_aggregate   <- NA_real_
    df$rmse_aggregate <- NA_real_
  } else {
    # Join aggregate values by subgroup_spec
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
#' \code{\link{summarytable}}. Results are grouped by
#' \code{subgroup_spec} and displayed with fixed-width columns for easy
#' reading. Aggregate values are printed once per subgroup block after the
#' per-stratum rows. The line width is computed dynamically from the content
#' so that separators always fit the output exactly.
#'
#' When \code{kable = TRUE}, the table is formatted using
#' \code{knitr::kable()}. For \code{kable_format = "html"}, additional
#' styling and row grouping are applied via \code{kableExtra} if available.
#'
#' @param x            An object of class \code{EvalMetrics}.
#' @param digits       Non-negative integer. Number of decimal places for
#'   NAUC, K-S distance, and RMSE values. Default \code{4}.
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
print.EvalMetrics <- function(x, digits = 4L, kable = FALSE,
                              kable_format = "pipe", ...) {

  kable_format <- match.arg(kable_format,
                            c("pipe", "simple", "latex", "html"))

  digits   <- as.integer(digits)
  tau_val  <- attr(x, "tau")
  time_var <- attr(x, "time_var")
  ev_var   <- attr(x, "event_var")
  n_meth   <- attr(x, "nauc_method_avg")
  k_meth   <- attr(x, "ks_method_avg")
  r_meth   <- attr(x, "rmse_method_avg")
  n_pts    <- attr(x, "n_points")

  # ---------------------------------------------------------------------------
  # Build header lines first so their widths can contribute to total_w
  # ---------------------------------------------------------------------------
  tau_str  <- if (is.null(tau_val)) "max observed" else
    formatC(tau_val, format = "f", digits = 1L)
  agg_str  <- sprintf("NAUC = %s | KS = %s | RMSE = %s",
                      n_meth, k_meth, r_meth)
  rmse_str <- sprintf("full observation range (#points = %d)", n_pts)

  # Fixed-label header lines (label width = 16 chars including trailing space)
  hdr_lines <- c(
    sprintf(" Time variable  : %s", time_var),
    sprintf(" Event variable : %s", ev_var),
    sprintf(" tau (NAUC/KS)  : %s", tau_str),
    sprintf(" RMSE range     : %s", rmse_str),
    sprintf(" Aggregation    : %s", agg_str)
  )

  # ---------------------------------------------------------------------------
  # Column widths computed from content
  # ---------------------------------------------------------------------------
  # Stratum label width: max of labels, "[Aggregate]", header text
  stratum_w <- max(
    nchar(x$stratum),
    nchar("[Aggregate]"),
    nchar("Subgroup / Stratum")
  ) + 1L

  # N columns: width 6 is sufficient for counts up to 999999
  n_w <- 6L

  # Metric column width: digits + 2 ("0.") but at least as wide as the header
  metric_val_w <- max(digits + 2L, nchar("NAUC"), nchar("KS"),
                      nchar("RMSE")) + 1L

  # Total width from the data columns:
  #  1 (leading space) + stratum_w + 2 (gap) + n_w + 1 (gap) + n_w +
  #  2 (gap) + metric_val_w + 1 (gap) + metric_val_w + 1 (gap) + metric_val_w
  col_total_w <- 1L + stratum_w + 2L + n_w + 1L + n_w +
    2L + metric_val_w + 1L + metric_val_w + 1L + metric_val_w

  # Ensure total width is never narrower than the widest header line
  total_w  <- max(col_total_w, max(nchar(hdr_lines)) + 1L)

  sep_eq   <- paste(rep("=", total_w), collapse = "")
  sep_dash <- paste(rep("-", total_w), collapse = "")

  # Format strings built from widths
  fmt_header <- sprintf(" %%-%ds %%%ds %%%ds  %%%ds %%%ds %%%ds\n",
                        stratum_w, n_w, n_w,
                        metric_val_w, metric_val_w, metric_val_w)
  fmt_row    <- sprintf("  %%-%ds %%%dd %%%dd  %%%d.%df %%%d.%df %%%d.%df\n",
                        stratum_w - 1L, n_w, n_w,
                        metric_val_w, digits,
                        metric_val_w, digits,
                        metric_val_w, digits)
  fmt_agg    <- sprintf("  %%-%ds %%%ds %%%ds  %%%d.%df %%%d.%df %%%d.%df\n",
                        stratum_w - 1L, n_w, n_w,
                        metric_val_w, digits,
                        metric_val_w, digits,
                        metric_val_w, digits)

  specs <- unique(x$subgroup_spec)

  # ---------------------------------------------------------------------------
  # Plain text output
  # ---------------------------------------------------------------------------
  if (!kable) {

    # Header block (plain text only)
    cat(sep_eq, "\n")
    cat(" EvalEmulateSurv: Emulation Evaluation Metrics\n")
    cat(sep_eq, "\n")
    for (hl in hdr_lines) cat(hl, "\n")
    cat(sep_dash, "\n")
    cat(sprintf(fmt_header,
                "Subgroup / Stratum", "N(g1)", "N(g2)",
                "NAUC", "KS", "RMSE"))
    cat(sep_dash, "\n")

    for (sp in specs) {
      rows <- x[x$subgroup_spec == sp, ]

      cat(sprintf(" [%s]\n", sp))

      for (i in seq_len(nrow(rows))) {
        r <- rows[i, ]
        cat(sprintf(fmt_row,
                    r$stratum,
                    r$n_group1, r$n_group2,
                    r$nauc, r$ks, r$rmse))
      }

      # Aggregate row (only when not Overall)
      if (sp != "Overall" && !is.na(rows$nauc_aggregate[1])) {
        cat(sprintf(fmt_agg,
                    "[Aggregate]", "", "",
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

    # Build a flat data.frame for kable, adding an [Aggregate] row per spec.
    # Metric values are pre-formatted as character strings via formatC() so
    # that knitr::kable() does not override the digit count with its own
    # automatic numeric formatting.
    fmt_metric <- function(v) formatC(v, format = "f", digits = digits)

    rows_list <- lapply(specs, function(sp) {
      rows <- x[x$subgroup_spec == sp, ]

      # Per-stratum rows
      df_str <- data.frame(
        Stratum  = rows$stratum,
        "N(g1)"  = rows$n_group1,
        "N(g2)"  = rows$n_group2,
        NAUC     = fmt_metric(rows$nauc),
        KS       = fmt_metric(rows$ks),
        RMSE     = fmt_metric(rows$rmse),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      # Append aggregate row if applicable
      if (sp != "Overall" && !is.na(rows$nauc_aggregate[1])) {
        agg_row <- data.frame(
          Stratum  = "[Aggregate]",
          "N(g1)"  = "",
          "N(g2)"  = "",
          NAUC     = fmt_metric(rows$nauc_aggregate[1]),
          KS       = fmt_metric(rows$ks_aggregate[1]),
          RMSE     = fmt_metric(rows$rmse_aggregate[1]),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        df_str <- rbind(df_str, agg_row)
      }

      df_str
    })

    kable_df           <- do.call(rbind, rows_list)
    rownames(kable_df) <- NULL

    # Build kable table; align metric columns right (character strings)
    tbl <- knitr::kable(
      kable_df,
      format  = kable_format,
      escape  = FALSE,
      align   = c("l", "r", "r", "r", "r", "r"),
      caption = sprintf(
        paste0("Emulation Evaluation Metrics | ",
               "time: %s | event: %s | tau: %s | ",
               "RMSE: full range (#pts=%d) | ",
               "Agg: NAUC=%s, KS=%s, RMSE=%s"),
        time_var, ev_var, tau_str, n_pts,
        n_meth, k_meth, r_meth
      )
    )

    # HTML: additional styling and subgroup row grouping via kableExtra
    if (kable_format == "html" &&
        requireNamespace("kableExtra", quietly = TRUE)) {

      tbl <- kableExtra::kable_styling(
        tbl,
        bootstrap_options = c("striped", "hover", "condensed", "responsive"),
        full_width        = FALSE,
        position          = "center"
      )

      # Add pack_rows() group labels per subgroup_spec
      row_offset <- 0L
      for (sp in specs) {
        rows    <- x[x$subgroup_spec == sp, ]
        n_rows  <- nrow(rows)
        has_agg <- sp != "Overall" && !is.na(rows$nauc_aggregate[1])
        n_block <- n_rows + if (has_agg) 1L else 0L

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
