#' Kaplan-Meier Plot with NAUC, K-S Distance, and RMSE Annotation
#'
#' @description
#' Draws Kaplan-Meier survival curves for two groups and annotates each panel
#' with NAUC, K-S distance, and RMSE. When \code{subgroup_var} is specified,
#' one panel is produced per stratum via \code{facet_wrap()}, and optionally
#' one for the overall population via \code{include_overall}. This function
#' supersedes the former \code{kmplot_by_subgroup()} function.
#'
#' @details
#' For each panel the following steps are performed:
#' \enumerate{
#'   \item KM step functions are estimated via \code{\link{km_est}} twice:
#'     once up to \code{tau} (for the ribbon and NAUC / K-S metrics) and
#'     once over the full observation range (for plotting).
#'   \item The area between the two KM curves from 0 to \code{tau} is shaded
#'     using exact rectangles derived from the step-function change points.
#'   \item A vertical dashed line marks the time of the maximum K-S gap
#'     within \eqn{[0, \tau]}.
#'   \item NAUC and K-S distance are computed over \eqn{[0, \tau]} via
#'     \code{\link{nauc}} and \code{\link{ksdist}}. RMSE is computed over the
#'     full observation range of each panel (\code{tau = NULL}) via
#'     \code{\link{rmse}}, following Lang et al. (2025).
#'   \item All three metrics are displayed in a monospace annotation block
#'     in the upper-right corner of each panel via \code{geom_text()}.
#' }
#'
#' All panels within a grid share the same x-axis range (the maximum observed
#' time across all strata), ensuring visual comparability.
#'
#' \strong{Subgroup specification via \code{subgroup_var}:}
#'
#' Each element of \code{subgroup_var} is either a single column name or an
#' interaction term using \code{":"} as separator (e.g., \code{"SEX:REGION"}).
#' For an interaction term, one panel is produced for each combination of the
#' factor levels (e.g., 2 x 3 = 6 panels for \code{"SEX:REGION"}).
#' N-way interactions are supported for any N.
#'
#' When multiple elements are supplied in \code{subgroup_var}, one
#' \code{ggplot} object is produced per element and returned as a named list.
#'
#' When \code{include_overall = TRUE}, an \code{"Overall"} panel is prepended
#' to each plot. When \code{subgroup_var = NULL}, \code{include_overall} is
#' ignored and a single overall panel is returned.
#'
#' \strong{Panel layout:}
#'
#' Use \code{ncol} to control the number of columns passed to
#' \code{facet_wrap()}. When \code{NULL}, \code{facet_wrap()} determines
#' the layout automatically.
#'
#' @param data_group1       Input for group 1 (e.g., original IPD). Accepted
#'   formats: a \code{data.frame} or tibble already in memory, a character
#'   path to a \code{.csv} file, or a character path to an \code{.RData} /
#'   \code{.rda} file containing a single object. Must contain columns
#'   specified by \code{time_var}, \code{event_var}, and all variables
#'   referenced in \code{subgroup_var}. See \code{\link{load_ipd}} for
#'   details.
#' @param data_group2       Input for group 2 (e.g., emulated IPD). Same
#'   accepted formats as \code{data_group1}.
#' @param time_var          Character. Name of the time-to-event column.
#'   Default \code{"TIME"}.
#' @param event_var         Character. Name of the event indicator column
#'   (1 = event, 0 = censored). Default \code{"EVENT"}.
#' @param tau               Numeric scalar. Truncation time used for the
#'   shaded ribbon, NAUC, and K-S distance. KM curves are plotted over the
#'   full observation range regardless of this value. RMSE is always computed
#'   over the full observation range of each panel. If \code{NULL} (default),
#'   the maximum observed time across both groups is used.
#' @param subgroup_var      Character vector or \code{NULL}. Each element is
#'   either a single column name (e.g., \code{"SEX"}) or an interaction term
#'   using \code{":"} as separator (e.g., \code{"SEX:REGION"}). N-way
#'   interactions are supported for any N. One \code{ggplot} object is
#'   produced per element and returned as a named list. If \code{NULL}
#'   (default), a single overall panel is produced.
#' @param include_overall   Logical. If \code{TRUE}, an \code{"Overall"} panel
#'   is prepended to each subgroup plot. Ignored when
#'   \code{subgroup_var = NULL}. Default \code{FALSE}.
#' @param ncol              Positive integer or \code{NULL}. Number of columns
#'   passed to \code{facet_wrap()}. Default \code{NULL} (automatic).
#' @param n_points          Positive integer. Number of evenly spaced time
#'   points used for RMSE calculation. Default \code{100L}.
#' @param label_group1      Character. Legend label for group 1.
#'   Default \code{"Original"}.
#' @param label_group2      Character. Legend label for group 2.
#'   Default \code{"Emulated"}.
#' @param color_group1      Character. Line colour for group 1.
#'   Default \code{"#004C97"} (dark blue).
#' @param color_group2      Character. Line colour for group 2.
#'   Default \code{"#F0B323"} (amber).
#' @param base_size         Numeric scalar. Base font size passed to
#'   \code{theme_bw()}. Default \code{14}.
#' @param metrics_text_size Numeric scalar. Font size (in ggplot2 units) for
#'   the annotation text block. Default \code{3}.
#'
#' @return A named list of \code{ggplot} objects. When
#'   \code{subgroup_var = NULL}, the list has a single element named
#'   \code{"Overall"}. When \code{subgroup_var} is specified, the list has
#'   one element per element of \code{subgroup_var}, named after each
#'   specification string. The list is returned invisibly.
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
#' @import ggplot2
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
#' plots <- kmplot(
#'   data_group1 = grp1_os,
#'   data_group2 = grp2_os,
#'   time_var    = "TIME",
#'   event_var   = "EVENT",
#'   tau         = 36
#' )
#' plots[["Overall"]]
#'
#' # -------------------------------------------------------------------
#' # 2. Single subgroup: REGION (3 panels), OS
#' # -------------------------------------------------------------------
#' plots_region <- kmplot(
#'   data_group1  = grp1_os,
#'   data_group2  = grp2_os,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 36,
#'   subgroup_var = "REGION",
#'   ncol         = 3
#' )
#' plots_region[["REGION"]]
#'
#' # -------------------------------------------------------------------
#' # 3. REGION with Overall prepended (4 panels)
#' # -------------------------------------------------------------------
#' plots_region_ov <- kmplot(
#'   data_group1     = grp1_os,
#'   data_group2     = grp2_os,
#'   time_var        = "TIME",
#'   event_var       = "EVENT",
#'   tau             = 36,
#'   subgroup_var    = "REGION",
#'   include_overall = TRUE,
#'   ncol            = 2
#' )
#' plots_region_ov[["REGION"]]
#'
#' # -------------------------------------------------------------------
#' # 4. Two-way interaction: SEX x REGION (6 panels), OS
#' # -------------------------------------------------------------------
#' plots_inter <- kmplot(
#'   data_group1  = grp1_os,
#'   data_group2  = grp2_os,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 36,
#'   subgroup_var = "SEX:REGION",
#'   ncol         = 3
#' )
#' plots_inter[["SEX:REGION"]]
#'
#' # -------------------------------------------------------------------
#' # 5. Multiple subgroup specs in one call
#' # -------------------------------------------------------------------
#' plots_multi <- kmplot(
#'   data_group1     = grp1_os,
#'   data_group2     = grp2_os,
#'   time_var        = "TIME",
#'   event_var       = "EVENT",
#'   tau             = 36,
#'   subgroup_var    = c("SEX", "REGION"),
#'   include_overall = TRUE,
#'   ncol            = 2
#' )
#' plots_multi[["SEX"]]
#' plots_multi[["REGION"]]
#'
#' # -------------------------------------------------------------------
#' # 6. PFS, SEX subgroup
#' # -------------------------------------------------------------------
#' plots_pfs <- kmplot(
#'   data_group1  = grp1_pfs,
#'   data_group2  = grp2_pfs,
#'   time_var     = "TIME",
#'   event_var    = "EVENT",
#'   tau          = 24,
#'   subgroup_var = "SEX",
#'   ncol         = 2
#' )
#' plots_pfs[["SEX"]]
#'
#' @export
kmplot <- function(data_group1,
                   data_group2,
                   time_var          = "TIME",
                   event_var         = "EVENT",
                   tau               = NULL,
                   subgroup_var      = NULL,
                   include_overall   = FALSE,
                   ncol              = NULL,
                   n_points          = 100L,
                   label_group1      = "Original",
                   label_group2      = "Emulated",
                   color_group1      = "#004C97",
                   color_group2      = "#F0B323",
                   base_size         = 14,
                   metrics_text_size = 3) {

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
  if (!is.logical(include_overall) || length(include_overall) != 1L)
    stop("'include_overall' must be a single logical value (TRUE or FALSE).")
  if (!is.numeric(n_points) || length(n_points) != 1L || n_points < 2L)
    stop("'n_points' must be a single integer >= 2.")
  if (!is.null(tau)) {
    if (!is.numeric(tau) || length(tau) != 1L || tau <= 0)
      stop("'tau' must be a single positive number.")
  }

  # ---------------------------------------------------------------------------
  # Internal helper: evaluate KM step function at time t
  # ---------------------------------------------------------------------------
  .km_eval <- function(km_df, t) {
    idx <- findInterval(t, km_df$time)
    if (idx == 0L) return(1.0)
    km_df$surv[idx]
  }

  # ---------------------------------------------------------------------------
  # Internal helper: build ribbon rectangles within [0, tau_val]
  # ---------------------------------------------------------------------------
  .build_ribbon_rects <- function(km1, km2, tau_val) {
    breaks <- sort(unique(c(km1$time, km2$time)))
    breaks <- breaks[breaks <= tau_val]
    breaks <- sort(unique(c(breaks, tau_val)))
    n      <- length(breaks)
    if (n < 2L) return(NULL)
    rects  <- vector("list", n - 1L)
    for (i in seq_len(n - 1L)) {
      s1 <- .km_eval(km1, breaks[i])
      s2 <- .km_eval(km2, breaks[i])
      rects[[i]] <- data.frame(
        xmin = breaks[i],   xmax = breaks[i + 1L],
        ymin = min(s1, s2), ymax = max(s1, s2)
      )
    }
    do.call(rbind, rects)
  }

  # ---------------------------------------------------------------------------
  # Internal helper: build annotation text string for one stratum
  # ---------------------------------------------------------------------------
  .make_annot_text <- function(tau_val, nauc_val, ks_val, t_max_gap,
                               rmse_val) {
    l_tau  <- "tau     "
    l_nauc <- "NAUC    "
    l_ks   <- "K-S dist"
    l_rmse <- "RMSE    "
    sep    <- ": "
    indent <- paste(rep(" ", nchar(l_ks) + nchar(sep)), collapse = "")
    paste(
      paste0(l_tau,  sep, formatC(tau_val,   format = "f", digits = 0,
                                  width = 6)),
      paste0(l_nauc, sep, formatC(nauc_val,  format = "f", digits = 4,
                                  width = 6)),
      paste0(l_ks,   sep, formatC(ks_val,    format = "f", digits = 4,
                                  width = 6)),
      paste0(indent, "(t = ",
             formatC(t_max_gap, format = "f", digits = 2), ")"),
      paste0(l_rmse, sep, formatC(rmse_val,  format = "f", digits = 4,
                                  width = 6)),
      paste0(indent, "(full range)"),
      sep = "\n"
    )
  }

  # ---------------------------------------------------------------------------
  # Internal helper: build a stratum character vector for one spec string.
  # Specs without ":" are single column names; specs with ":" are N-way
  # interaction terms. Stratum labels use " / " as separator.
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
  # Internal helper: build a single ggplot for one specification.
  #
  # All strata data are combined into one data.frame and facet_wrap() is used
  # to produce one panel per stratum. This avoids patchwork entirely, so no
  # guide-collection warnings can arise. The shared legend is produced
  # naturally by ggplot2 and positioned at the bottom via theme().
  # ---------------------------------------------------------------------------
  .build_plot <- function(spec) {

    # Determine stratum labels for each observation
    if (spec == "Overall") {
      strata_d1 <- rep("Overall", nrow(data_group1))
      strata_d2 <- rep("Overall", nrow(data_group2))
    } else {
      strata_d1 <- .make_strata_col(data_group1, spec)
      strata_d2 <- .make_strata_col(data_group2, spec)
    }

    strata <- sort(unique(c(strata_d1, strata_d2)))

    # Prepend Overall stratum when include_overall = TRUE and spec != "Overall"
    if (spec != "Overall" && include_overall) {
      strata <- c("Overall", strata)
    }

    # Global x-axis upper limit: max t_max_obs across all strata in this plot
    x_max_global <- max(
      c(data_group1[[time_var]], data_group2[[time_var]]),
      na.rm = TRUE
    )

    # Build per-stratum data for step curves, ribbons, and annotations
    step_list   <- list()
    ribbon_list <- list()
    annot_list  <- list()

    for (s in strata) {
      if (s == "Overall") {
        d1 <- data_group1
        d2 <- data_group2
      } else {
        d1 <- data_group1[strata_d1 == s, ]
        d2 <- data_group2[strata_d2 == s, ]
      }

      if (nrow(d1) == 0L || nrow(d2) == 0L) {
        warning(sprintf(
          "Spec '%s', stratum '%s': empty subset in one group. Panel skipped.",
          spec, s
        ))
        next
      }

      t_max_obs <- max(c(d1[[time_var]], d2[[time_var]]), na.rm = TRUE)
      tau_s     <- if (is.null(tau)) t_max_obs else tau

      # KM for ribbon / metrics (restricted to tau_s)
      km_tau  <- km_est(d1, d2, time_var = time_var,
                        event_var = event_var, tau = tau_s)
      # KM for plotting (full range)
      km_full <- km_est(d1, d2, time_var = time_var,
                        event_var = event_var, tau = t_max_obs)

      # Metrics
      res_nauc <- nauc(d1, d2, time_var = time_var,
                       event_var = event_var, tau = tau_s)
      res_ks   <- ksdist(d1, d2, time_var = time_var,
                         event_var = event_var, tau = tau_s)
      res_rmse <- rmse(d1, d2, time_var = time_var,
                       event_var = event_var, tau = NULL,
                       n_points = as.integer(n_points))

      nauc_val  <- res_nauc$aggregate$nauc_aggregate[1]
      ks_val    <- res_ks$aggregate$ks_aggregate[1]
      t_max_gap <- res_ks$by_stratum$t_max_gap[1]
      rmse_val  <- res_rmse$aggregate$rmse_aggregate[1]

      # Step curve data (full range)
      step_list[[s]] <- rbind(
        data.frame(time    = km_full$km1$time,
                   surv    = km_full$km1$surv,
                   group   = label_group1,
                   stratum = s,
                   stringsAsFactors = FALSE),
        data.frame(time    = km_full$km2$time,
                   surv    = km_full$km2$surv,
                   group   = label_group2,
                   stratum = s,
                   stringsAsFactors = FALSE)
      )

      # Ribbon data
      rects <- .build_ribbon_rects(km_tau$km1, km_tau$km2, tau_s)
      if (!is.null(rects)) {
        rects$stratum    <- s
        ribbon_list[[s]] <- rects
      }

      # Vertical line and annotation data
      annot_list[[s]] <- data.frame(
        stratum   = s,
        t_max_gap = t_max_gap,
        x_annot   = x_max_global,
        label     = .make_annot_text(tau_s, nauc_val, ks_val,
                                     t_max_gap, rmse_val),
        stringsAsFactors = FALSE
      )
    }

    if (length(step_list) == 0L)
      stop(sprintf("No valid panels could be built for spec '%s'.", spec))

    # Combine all strata into single data.frames
    step_df   <- do.call(rbind, step_list)
    ribbon_df <- do.call(rbind, ribbon_list)
    annot_df  <- do.call(rbind, annot_list)

    # Fix factor order so that panels appear in the correct sequence
    step_df$stratum   <- factor(step_df$stratum,   levels = strata)
    ribbon_df$stratum <- factor(ribbon_df$stratum, levels = strata)
    annot_df$stratum  <- factor(annot_df$stratum,  levels = strata)

    step_df$group <- factor(step_df$group,
                            levels = c(label_group1, label_group2))

    # Build ggplot with facet_wrap
    ggplot2::ggplot() +
      ggplot2::geom_rect(
        data    = ribbon_df,
        mapping = ggplot2::aes(
          xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
        ),
        fill   = "grey70",
        alpha  = 0.6,
        colour = NA,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_step(
        data      = step_df,
        mapping   = ggplot2::aes(x = time, y = surv, colour = group),
        linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        data     = annot_df,
        mapping  = ggplot2::aes(xintercept = t_max_gap),
        linetype = "dashed",
        colour   = "grey30",
        linewidth = 0.6,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data    = annot_df,
        mapping = ggplot2::aes(x = x_annot, y = 0.98, label = label),
        hjust   = 1,
        vjust   = 1,
        size    = metrics_text_size,
        family  = "mono",
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = setNames(c(color_group1, color_group2),
                          c(label_group1, label_group2))
      ) +
      ggplot2::scale_x_continuous(limits = c(0, x_max_global)) +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, by = 0.2)) +
      ggplot2::facet_wrap(~ stratum, ncol = ncol) +
      ggplot2::labs(x = "Time", y = "Survival Probability") +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        legend.position  = "bottom",
        legend.text      = ggplot2::element_text(size = base_size * 0.8),
        legend.title     = ggplot2::element_blank(),
        legend.key.width = ggplot2::unit(2, "cm"),
        panel.grid.minor = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(size = base_size * 0.9)
      )
  }

  # ---------------------------------------------------------------------------
  # Step 3: Build list of specifications to process
  # ---------------------------------------------------------------------------
  specs <- if (is.null(subgroup_var)) "Overall" else subgroup_var

  # ---------------------------------------------------------------------------
  # Step 4: Build one ggplot per specification and return as named list
  # ---------------------------------------------------------------------------
  plot_list        <- lapply(specs, .build_plot)
  names(plot_list) <- specs

  invisible(plot_list)
}
