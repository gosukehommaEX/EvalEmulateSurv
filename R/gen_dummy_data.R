#' Generate Dummy IPD for Package Examples and Testing
#'
#' @description
#' Generates a dummy individual patient data (IPD) data frame containing
#' survival endpoints (OS and PFS), treatment group, and covariates (SEX
#' and REGION). Survival times are drawn from exponential distributions.
#' The function is intended to provide a self-contained example dataset
#' for package documentation and unit testing.
#'
#' @details
#' For each combination of endpoint and treatment group, survival times are
#' generated independently from an exponential distribution with the
#' corresponding rate parameter. Event indicators are drawn from a Bernoulli
#' distribution with probability \code{event_prob}. Covariate values (SEX,
#' REGION) are assigned by simple random sampling with replacement using the
#' proportions specified in \code{sex_prob} and \code{region_prob}.
#'
#' The returned data frame has one row per patient per endpoint. The
#' \code{TYPE} column distinguishes the two endpoints (\code{"OS"} and
#' \code{"PFS"}), enabling downstream functions (\code{km_est()},
#' \code{nauc()}, \code{ks_distance()}, \code{rmse()}) to be applied
#' after filtering by \code{TYPE}.
#'
#' @param n_drug     Integer. Number of patients in the Drug group.
#'   Default \code{80}.
#' @param n_control  Integer. Number of patients in the Control group.
#'   Default \code{80}.
#' @param rate_os_drug     Numeric. Exponential rate for OS in the Drug
#'   group. Default \code{0.030}.
#' @param rate_os_control  Numeric. Exponential rate for OS in the Control
#'   group. Default \code{0.042}.
#' @param rate_pfs_drug    Numeric. Exponential rate for PFS in the Drug
#'   group. Default \code{0.055}.
#' @param rate_pfs_control Numeric. Exponential rate for PFS in the Control
#'   group. Default \code{0.075}.
#' @param event_prob Numeric in \code{[0, 1]}. Probability that an event
#'   (not censoring) occurred. Applied identically across all groups and
#'   endpoints. Default \code{0.85}.
#' @param sex_prob   Named numeric vector of length 2 giving sampling
#'   probabilities for \code{c("M", "F")}. Values need not sum to 1 (they
#'   are normalised internally). Default \code{c(M = 0.5, F = 0.5)}.
#' @param region_prob Named numeric vector of length 3 giving sampling
#'   probabilities for \code{c("US", "EU", "OTHERS")}. Values need not
#'   sum to 1. Default \code{c(US = 0.4, EU = 0.35, OTHERS = 0.25)}.
#' @param seed Integer or \code{NULL}. Random seed passed to
#'   \code{set.seed()} for reproducibility. Set to \code{NULL} to disable.
#'   Default \code{42}.
#'
#' @return A data frame with \code{(n_drug + n_control) * 2} rows and the
#'   following columns:
#' \describe{
#'   \item{\code{ID}}{Integer patient identifier (unique within each
#'     \code{GROUP}).}
#'   \item{\code{GROUP}}{Character. Treatment group: \code{"Drug"} or
#'     \code{"Control"}.}
#'   \item{\code{TYPE}}{Character. Endpoint name: \code{"OS"} or
#'     \code{"PFS"}.}
#'   \item{\code{TIME}}{Numeric. Observed time (survival or censoring
#'     time), drawn from an exponential distribution.}
#'   \item{\code{EVENT}}{Integer. Event indicator: \code{1} = event
#'     occurred, \code{0} = censored.}
#'   \item{\code{SEX}}{Character. Patient sex: \code{"M"} or \code{"F"}.}
#'   \item{\code{REGION}}{Character. Geographic region: \code{"US"},
#'     \code{"EU"}, or \code{"OTHERS"}.}
#' }
#'
#' @import stats
#'
#' @examples
#' # Default settings
#' ipd <- gen_dummy_data()
#' head(ipd)
#' table(ipd$GROUP, ipd$TYPE)
#'
#' # Filter to OS endpoint for km_est()
#' ipd_os <- ipd[ipd$TYPE == "OS", ]
#' orig    <- ipd_os[ipd_os$GROUP == "Drug", ]
#' emu     <- ipd_os[ipd_os$GROUP == "Control", ]
#'
#' # Smaller dataset with custom seed
#' ipd_small <- gen_dummy_data(n_drug = 40, n_control = 40, seed = 123)
#' nrow(ipd_small)
#'
#' @export
gen_dummy_data <- function(n_drug          = 80L,
                           n_control       = 80L,
                           rate_os_drug    = 0.030,
                           rate_os_control = 0.042,
                           rate_pfs_drug   = 0.055,
                           rate_pfs_control = 0.075,
                           event_prob      = 0.85,
                           sex_prob        = c(M = 0.5, F = 0.5),
                           region_prob     = c(US = 0.4, EU = 0.35,
                                               OTHERS = 0.25),
                           seed            = 42L) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.numeric(n_drug) || length(n_drug) != 1L || n_drug < 1L)
    stop("'n_drug' must be a single positive integer.")
  if (!is.numeric(n_control) || length(n_control) != 1L || n_control < 1L)
    stop("'n_control' must be a single positive integer.")

  for (rate_name in c("rate_os_drug", "rate_os_control",
                      "rate_pfs_drug", "rate_pfs_control")) {
    val <- get(rate_name)
    if (!is.numeric(val) || length(val) != 1L || val <= 0)
      stop(sprintf("'%s' must be a single positive number.", rate_name))
  }

  if (!is.numeric(event_prob) || length(event_prob) != 1L ||
      event_prob < 0 || event_prob > 1)
    stop("'event_prob' must be a single number in [0, 1].")

  if (!is.numeric(sex_prob) || length(sex_prob) != 2L || any(sex_prob < 0))
    stop("'sex_prob' must be a non-negative numeric vector of length 2.")
  if (!is.numeric(region_prob) || length(region_prob) != 3L ||
      any(region_prob < 0))
    stop("'region_prob' must be a non-negative numeric vector of length 3.")

  # ---------------------------------------------------------------------------
  # Set random seed
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  n_drug    <- as.integer(n_drug)
  n_control <- as.integer(n_control)

  # ---------------------------------------------------------------------------
  # Internal helper: build one-endpoint block for a single group
  # ---------------------------------------------------------------------------
  .make_block <- function(n, rate, group_label, type_label) {
    data.frame(
      ID    = seq_len(n),
      GROUP = group_label,
      TYPE  = type_label,
      TIME  = rexp(n, rate = rate),
      EVENT = rbinom(n, 1L, event_prob),
      stringsAsFactors = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Generate survival data for each group x endpoint combination
  # ---------------------------------------------------------------------------
  drug_os    <- .make_block(n_drug,    rate_os_drug,     "Drug",    "OS")
  control_os <- .make_block(n_control, rate_os_control,  "Control", "OS")
  drug_pfs   <- .make_block(n_drug,    rate_pfs_drug,    "Drug",    "PFS")
  control_pfs <- .make_block(n_control, rate_pfs_control, "Control", "PFS")

  ipd <- rbind(drug_os, control_os, drug_pfs, control_pfs)

  # ---------------------------------------------------------------------------
  # Assign covariates (SEX and REGION) consistent within patient x group
  # Covariates are assigned per patient (ID x GROUP combination) and then
  # merged back so that OS and PFS rows for the same patient share the same
  # covariate values.
  # ---------------------------------------------------------------------------
  sex_levels    <- c("M", "F")
  region_levels <- c("US", "EU", "OTHERS")

  # Unique patients per group
  drug_patients    <- seq_len(n_drug)
  control_patients <- seq_len(n_control)
  n_all            <- n_drug + n_control

  sex_all    <- sample(sex_levels,    n_all, replace = TRUE,
                       prob = sex_prob / sum(sex_prob))
  region_all <- sample(region_levels, n_all, replace = TRUE,
                       prob = region_prob / sum(region_prob))

  covar_drug <- data.frame(
    ID     = drug_patients,
    GROUP  = "Drug",
    SEX    = sex_all[seq_len(n_drug)],
    REGION = region_all[seq_len(n_drug)],
    stringsAsFactors = FALSE
  )
  covar_control <- data.frame(
    ID     = control_patients,
    GROUP  = "Control",
    SEX    = sex_all[n_drug + seq_len(n_control)],
    REGION = region_all[n_drug + seq_len(n_control)],
    stringsAsFactors = FALSE
  )
  covar <- rbind(covar_drug, covar_control)

  # ---------------------------------------------------------------------------
  # Merge covariates and return
  # ---------------------------------------------------------------------------
  ipd <- merge(ipd, covar, by = c("ID", "GROUP"), all.x = TRUE)

  # Restore row order: Drug OS, Control OS, Drug PFS, Control PFS
  ipd <- ipd[order(
    match(ipd$TYPE,  c("OS", "PFS")),
    match(ipd$GROUP, c("Drug", "Control")),
    ipd$ID
  ), ]
  row.names(ipd) <- NULL

  ipd
}
