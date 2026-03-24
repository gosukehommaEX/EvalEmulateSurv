#' Load Individual Patient Data from Various Sources
#'
#' @description
#' A utility function that accepts individual patient data (IPD) from multiple
#' input formats and returns a plain \code{data.frame}. Supported formats are:
#' \itemize{
#'   \item A \code{data.frame} or \code{tbl_df} (tibble) already in memory.
#'   \item A path to a \code{.csv} file (read with \code{read.csv()}).
#'   \item A path to an \code{.RData} / \code{.rda} file containing exactly
#'     one object, or the name of the target object within a multi-object
#'     \code{.RData} file supplied via \code{rdata_object}.
#' }
#' The function is intended as a pre-processing helper for
#' \code{nauc()}, \code{ks_distance()}, \code{rmse()}, and \code{km_est()},
#' which all expect a \code{data.frame} for their \code{data_group1} and
#' \code{data_group2} arguments.
#'
#' @details
#' \strong{Dispatch logic:}
#' \enumerate{
#'   \item If \code{x} is already a \code{data.frame} or tibble, it is
#'     coerced to a plain \code{data.frame} and returned immediately.
#'   \item If \code{x} is a character string ending in \code{.csv}
#'     (case-insensitive), the file is read with \code{read.csv()}.
#'     Additional arguments can be passed via \code{...} (e.g.,
#'     \code{stringsAsFactors = FALSE}).
#'   \item If \code{x} is a character string ending in \code{.RData} or
#'     \code{.rda} (case-insensitive), the file is loaded into an isolated
#'     temporary environment. When \code{rdata_object = NULL} the file must
#'     contain exactly one object; otherwise \code{rdata_object} names the
#'     target object.
#'   \item All other inputs raise an informative error.
#' }
#'
#' @param x An object to load. Accepted types:
#'   \itemize{
#'     \item \code{data.frame} or \code{tbl_df}: returned as-is (coerced to
#'       \code{data.frame}).
#'     \item \code{character(1)}: path to a \code{.csv}, \code{.RData}, or
#'       \code{.rda} file.
#'   }
#' @param rdata_object \code{NULL} (default) or a character string giving the
#'   name of the object to extract from an \code{.RData} / \code{.rda} file
#'   that contains more than one object. Ignored for \code{.csv} input and
#'   in-memory objects.
#' @param ... Additional arguments forwarded to \code{read.csv()} when
#'   \code{x} is a path to a \code{.csv} file. Ignored for other input types.
#'
#' @return A \code{data.frame}.
#'
#' @seealso \code{\link{km_est}}, \code{\link{nauc}},
#'   \code{\link{ks_distance}}, \code{\link{rmse}}
#'
#' @import utils
#'
#' @examples
#' # -------------------------------------------------------------------
#' # 1. In-memory data.frame (most common usage)
#' # -------------------------------------------------------------------
#' ipd <- gen_dummy_data()
#' orig <- load_ipd(ipd[ipd$GROUP == "Drug" & ipd$TYPE == "OS", ])
#' class(orig)  # "data.frame"
#'
#' # -------------------------------------------------------------------
#' # 2. In-memory tibble
#' # -------------------------------------------------------------------
#' if (requireNamespace("dplyr", quietly = TRUE)) {
#'   library(dplyr)
#'   ipd_tbl <- dplyr::as_tibble(gen_dummy_data())
#'   orig_tbl <- load_ipd(ipd_tbl[ipd_tbl$GROUP == "Drug" &
#'                                  ipd_tbl$TYPE == "OS", ])
#'   class(orig_tbl)  # "data.frame"
#' }
#'
#' # -------------------------------------------------------------------
#' # 3. CSV file
#' # -------------------------------------------------------------------
#' tmp_csv <- tempfile(fileext = ".csv")
#' write.csv(gen_dummy_data(), tmp_csv, row.names = FALSE)
#' ipd_csv <- load_ipd(tmp_csv)
#' head(ipd_csv)
#' unlink(tmp_csv)
#'
#' # -------------------------------------------------------------------
#' # 4. RData file containing a single object
#' # -------------------------------------------------------------------
#' tmp_rda <- tempfile(fileext = ".RData")
#' ipd_orig <- gen_dummy_data()
#' save(ipd_orig, file = tmp_rda)
#' ipd_rda <- load_ipd(tmp_rda)
#' head(ipd_rda)
#' unlink(tmp_rda)
#'
#' # -------------------------------------------------------------------
#' # 5. RData file containing multiple objects (rdata_object required)
#' # -------------------------------------------------------------------
#' tmp_rda2 <- tempfile(fileext = ".RData")
#' ipd_drug    <- gen_dummy_data()[gen_dummy_data()$GROUP == "Drug", ]
#' ipd_control <- gen_dummy_data()[gen_dummy_data()$GROUP == "Control", ]
#' save(ipd_drug, ipd_control, file = tmp_rda2)
#' orig_multi <- load_ipd(tmp_rda2, rdata_object = "ipd_drug")
#' head(orig_multi)
#' unlink(tmp_rda2)
#'
#' @export
load_ipd <- function(x, rdata_object = NULL, ...) {

  # ---------------------------------------------------------------------------
  # Branch 1: already a data.frame or tibble
  # ---------------------------------------------------------------------------
  if (is.data.frame(x)) {
    return(as.data.frame(x))
  }

  # ---------------------------------------------------------------------------
  # Branch 2: file path (character scalar)
  # ---------------------------------------------------------------------------
  if (is.character(x) && length(x) == 1L) {

    if (!file.exists(x))
      stop(sprintf("File not found: '%s'", x))

    ext <- tolower(tools::file_ext(x))

    # -- CSV ------------------------------------------------------------------
    if (ext == "csv") {
      out <- utils::read.csv(x, ...)
      if (!is.data.frame(out))
        stop("'read.csv()' did not return a data.frame.")
      return(out)
    }

    # -- RData / rda ----------------------------------------------------------
    if (ext %in% c("rdata", "rda")) {
      env   <- new.env(parent = emptyenv())
      nms   <- load(x, envir = env)

      if (is.null(rdata_object)) {
        # Expect exactly one object
        if (length(nms) != 1L)
          stop(sprintf(
            paste0("'%s' contains %d objects (%s). ",
                   "Specify the target via 'rdata_object'."),
            x, length(nms), paste(nms, collapse = ", ")
          ))
        obj <- get(nms[1L], envir = env)
      } else {
        # User-specified object name
        if (!rdata_object %in% nms)
          stop(sprintf(
            "'rdata_object = \"%s\"' not found in '%s'. Available: %s",
            rdata_object, x, paste(nms, collapse = ", ")
          ))
        obj <- get(rdata_object, envir = env)
      }

      if (!is.data.frame(obj))
        stop(sprintf(
          "The loaded object is not a data.frame (class: %s).",
          paste(class(obj), collapse = ", ")
        ))
      return(as.data.frame(obj))
    }

    # -- Unsupported extension ------------------------------------------------
    stop(sprintf(
      paste0("Unsupported file extension '.%s'. ",
             "Accepted extensions: .csv, .RData, .rda."),
      ext
    ))
  }

  # ---------------------------------------------------------------------------
  # Branch 3: unsupported type
  # ---------------------------------------------------------------------------
  stop(sprintf(
    paste0("'x' must be a data.frame, tibble, or a character path to a ",
           ".csv / .RData / .rda file (got class: %s)."),
    paste(class(x), collapse = ", ")
  ))
}
